import torch
import torch.nn as nn
import torch.autograd as autograd
import torch.optim as optim
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams

import math

import torch.optim.adamw
from torch.utils.data import TensorDataset, random_split, DataLoader
from scipy.integrate import solve_ivp
from scipy.special import comb

torch.set_default_dtype(torch.float32)
torch.backends.mps.is_available()

rcParams['pdf.fonttype'] = 42  # Ensures TrueType (Type 1) fonts in PDFs
rcParams['ps.fonttype'] = 42
rcParams['font.family'] = 'Times New Roman'


sms_del = 1e-6
sms_cm  = 1.0
sms_or  = 1

def cb(n, k):
    result = 1
    for i in range(k):
        result *= (n - i) / (i + 1)
    return result

def sig(r):
    return math.sqrt(r**2 + sms_del**2) - sms_del

def smstep(sigma):
    b = sig(sms_cm)
    temp = 0
    if sigma <= b:
        for k in range(sms_or + 1):
            temp += comb(sms_or + k, k) * comb(2 * sms_or + 1, sms_or - k) * (-sigma / b) ** k
        temp *= (sigma / b) ** (sms_or + 1)
    else:
        temp = 1.0
    return temp

def d_smstep(sigma):
    b = sig(sms_cm)
    temp = 0
    if sigma <= b:
        for k in range(sms_or + 1):
            temp += cb(sms_or + k, k) * cb(2 * sms_or + 1, sms_or - k) * (-1) ** k * (sigma / b) ** (sms_or + k) * (sms_or + k + 1) / b
    return temp

# ==================================
# Define model 
# ==================================
class ReOrPHNN(nn.Module):
    """
    Attributes:
    - `H_net`: Neural network parameterizing H(x), outputs a scalar.
    - `V`: Transformation matrix.
    """
    def __init__(self, H_net, V):
        super(ReOrPHNN, self).__init__()
        self.H_net = H_net  # Neural network parameterizing H(x)
        self.V = V          # Constant input matrix
        self.gamma = torch.nn.Parameter(torch.tensor(1.0, requires_grad=True))
    
    def forward(self, x, u):
        n, r = self.V.shape
        N    = int(n/2)

        # Construct the J matrix (n x n)
        J = torch.cat([
            torch.cat([torch.zeros((N, N)), torch.eye(N)], dim=1),
            torch.cat([-torch.eye(N), torch.zeros((N, N))], dim=1)
        ], dim=0)
        Jr = self.V.T @ J @ self.V

        # Define the R matrix (n x n) with block diagonal structure
        self.gamma.data.clamp_(0.0, 1.0)

        R = torch.cat([
            torch.cat([torch.zeros((N, N)), torch.zeros((N, N))], dim=1),
            torch.cat([torch.zeros((N, N)), self.gamma*torch.eye(N)],  dim=1)], dim=0)
        Rr = self.V.T @ R @ self.V

        # Define the B matrix (n x 1)
        e1 = torch.zeros(N)
        e1[0] = 1  # e1 is the first standard basis vector
        B = torch.cat([torch.zeros(N), e1]).view(n, 1)
        Br = self.V.T @ B

        # Ensure x requires gradients for computing dH/dx
        x   = x.clone().detach().requires_grad_(True)  # Ensure x has gradients

        # Compute dH/dx using autograd
        mag  = torch.norm(x, p=2)
        NN_x = self.H_net(x)
        # H_x  = NN_x*smstep(sig(mag))
        # dHdx = autograd.grad(NN_x, x, create_graph=True)[0]*smstep(sig(mag)) + x*NN_x*d_smstep(sig(mag))/math.sqrt(mag**2 + sms_del**2)

        H_x  = sig(NN_x)*smstep(sig(mag))
        dHdx = autograd.grad(NN_x, x, create_graph=True)[0]*smstep(sig(mag))*NN_x/math.sqrt(NN_x**2 + sms_del**2) + x*NN_x*d_smstep(sig(mag))/math.sqrt(mag**2 + sms_del**2)

        dxdt = (Jr - Rr) @ dHdx + Br @ u  # Compute dxdt
        y    = (Br.T @ dHdx).squeeze()
        return dxdt, y, H_x
    
    def forward_batch(self, x, u):
        """
        Computes the forward pass for batched PH system.
        Args:
        - `x`: State matrix, where each column is a state vector.
        - `u`: Input vector, where each column is an scalar.
        Returns:
        - The time derivative of the state matrix, x'.
        """
        batch_size = x.shape[0]
        return torch.stack([self.forward(x[i, :], u[i, :])[0] for i in range(batch_size)], dim=0)

# ==================================
# Load the training data 
# ==================================
data_raw = np.load(os.path.join("..","data", "TodaLat_data_train.npz"))
DX    = data_raw["Xs"]
DY    = data_raw["Ys"]
DU    = data_raw["Us"]
DXdot = data_raw["Xdots"]
V     = data_raw["Vtrans"]
splt  = data_raw["Ts"]

DX    = V.T @ DX
DXdot = V.T @ DXdot

# Convert numpy arrays to torch tensors
DX    = torch.tensor(DX.T,              dtype = torch.float32)   # Transpose DX
DU    = torch.tensor(DU.reshape(-1, 1), dtype = torch.float32)   # Transpose DU
DY    = torch.tensor(DY.reshape(-1, 1), dtype = torch.float32)   # Transpose DU
DXdot = torch.tensor(DXdot.T,           dtype = torch.float32)   # Transpose DXdot
splt  = torch.tensor(splt, dtype = torch.float32)
V     = torch.tensor(V,    dtype = torch.float32)

# Create a dataset
dataset = TensorDataset(DX, DU, DXdot)

# Define train-validation split ratio
val_ratio  = 0.2
val_size   = int(len(dataset) * val_ratio)
train_size = len(dataset) - val_size

# Perform random split
train_dataset, val_dataset = random_split(dataset, [train_size, val_size])

# Create DataLoaders
batch_size   = 8
train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
val_loader   = DataLoader(val_dataset,   batch_size=batch_size, shuffle=False)

# ==================================
# Construct the network architecture
# ==================================
nu = len(DU.shape)
nx = DX.shape[1]

# Define H_net using nn.Sequential
nncells = 20
H_net = nn.Sequential(
    nn.Linear(nx, nncells),      # Input layer: nx -> nncells
    nn.Tanh(),                   # Activation function
    nn.Linear(nncells, nncells), # Hidden layer: nncells -> nncells
    nn.Tanh(),                   # Activation function
    nn.Linear(nncells, 1),       # Output layer: nncells -> 1 (scalar output)
    # nn.Softmax(dim=0)
)   
mor_PH = ReOrPHNN(H_net, V)


# ==================================
# TRAINING
# ==================================

# Define training parameters
num_epochs      = 5000  # Max number of epochs
learning_rate   = 2e-3
patience        = 200   # Number of epochs to wait before stopping
min_delta       = 1e-7  # Minimum improvement to reset patience
scheduler_step  = 500   # Reduce LR every 500 epochs
scheduler_gamma = 0.95  # LR reduction factor

# Define optimizer and loss function
# torch.optim.SGD
optimizer = torch.optim.AdamW(mor_PH.parameters(), lr=learning_rate)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=scheduler_step, gamma=scheduler_gamma)
criterion = nn.MSELoss()

# Initialize lists to store losses
train_losses = []
val_losses   = []

# Early stopping variables
best_val_loss = float("inf")
patience_counter = 0

print("\nStarting training...\n")

# Training loop
for epoch in range(1, num_epochs + 1):
    mor_PH.train()  # Set model to training mode
    total_train_loss = 0.0

    # Training phase
    for x_batch, u_batch, x_dot_batch in train_loader:
        optimizer.zero_grad()  # Reset gradients
        
        # Compute predicted x_dot
        x_dot_pred = mor_PH.forward_batch(x_batch, u_batch)

        # Compute loss
        loss = criterion(x_dot_pred, x_dot_batch)

        # Constraint loss
        total_loss = loss   # Combine with main loss
        total_loss.backward()

        optimizer.step()  # Update model parameters

        total_train_loss += loss.item()

    # Average training loss
    avg_train_loss = total_train_loss / len(train_loader)
    train_losses.append(avg_train_loss)

    # Validation phase
    mor_PH.eval()  # Set model to evaluation mode
    total_val_loss = 0.0

    for x_val, u_val, x_dot_val in val_loader:
        x_dot_pred_val = mor_PH.forward_batch(x_val, u_val)  # Compute predicted x_dot
        val_loss = criterion(x_dot_pred_val, x_dot_val)
        total_val_loss += val_loss.item()

    # Average validation loss
    avg_val_loss = total_val_loss / len(val_loader)
    val_losses.append(avg_val_loss)

    # Learning rate scheduler step
    scheduler.step()

    # Early stopping check
    if avg_val_loss < best_val_loss - min_delta:
        best_val_loss = avg_val_loss
        patience_counter = 0  # Reset patience counter
    else:
        patience_counter += 1

    if patience_counter >= patience:
        print(f"\nEarly stopping triggered at epoch {epoch} (patience {patience} reached)\n")
        break  # Stop training

    # Print progress summary every 50 epochs
    if epoch % 50 == 0:
        print(f"Epoch {epoch}: Train Loss = {avg_train_loss:.6f}, Val Loss = {avg_val_loss:.6f}, LR = {scheduler.get_last_lr()[0]:.6e}")

print("\nTraining complete!\n")

# Plot training & validation losses
plt.figure(figsize=(8, 5))
plt.plot(train_losses, label="Training Loss", color="blue")
plt.plot(val_losses, label="Validation Loss", color="red")
plt.xlabel("Epochs")
plt.ylabel("Loss")
plt.title("Training and Validation Loss Over Epochs")
plt.legend()
plt.grid(True)
plt.show()




# ==================================
# Load the test data
# ==================================
data_raw = np.load(os.path.join("..","data", "TodaLat_data_test.npz"))
test_DX, test_DU, test_DXdot = data_raw["Xs"], data_raw["Us"], data_raw["Xdots"]
test_splt = data_raw["Ts"]
test_DY   = data_raw["Ys"]
test_H    = data_raw["Hs"]

# Convert numpy arrays to torch tensors with necessary transformations
test_DX    = V.T.numpy() @ test_DX
test_DXdot = V.T.numpy() @ test_DXdot

test_DX    = torch.tensor(test_DX.T,    dtype=torch.float32)      
test_DXdot = torch.tensor(test_DXdot.T, dtype=torch.float32)
test_DU    = torch.tensor(test_DU.reshape(-1, 1), dtype=torch.float32)

# Create a dataset for evaluation
test_dataset = TensorDataset(test_DX, test_DU, test_DXdot)
# Create DataLoader for batch evaluation
test_loader  = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

# ================================
# Model Evaluation on Test Dataset
# ================================
mor_PH.eval()  # Ensure model is in evaluation mode

test_loss = 0.0

for x_test, u_test, x_dot_test in test_loader:
    x_dot_pred_test = mor_PH.forward_batch(x_test, u_test)  
    loss = criterion(x_dot_pred_test, x_dot_test)
    test_loss += loss.item()

test_loss /= len(test_loader)

print(f"\nTest MSE Loss: {test_loss:.6f}")

# ================================
# Model Evaluation by solving ODE
# ================================

# Extract initial state x(0) and input sequence U from the test dataset
x0      = test_DX[0].numpy()
U_seq   = test_DU.numpy()     # Input sequence
T_final = test_splt[-1]      # Final time (seconds)
dt      = test_splt[1] - test_splt[0]  # Time step

# Time vector for interpolation
t_eval = np.arange(0, T_final, dt)

# Define ODE system using the trained model
def ode_system(t, x):
    # Find the closest input value based on time index
    u_index = min(int(t / dt), len(U_seq) - 1)
    u = torch.tensor(U_seq[u_index], dtype=torch.float32).reshape(-1, 1)

    # Convert x to tensor
    x_tensor = torch.tensor(x, dtype=torch.float32).reshape(-1, 1)

    # Compute x' using the trained model (without torch.no_grad())
    x_dot = mor_PH.forward(x_tensor.squeeze(1), u.squeeze(1))[0]

    # Convert tensor to numpy array
    return x_dot.detach().numpy().flatten()  

# Solve the ODE using solve_ivp

sol = solve_ivp(ode_system, [test_splt[0], T_final], x0, method='RK45', t_eval=t_eval)

# Extract simulation results
simulated_X   = sol.y.T
true_X        = data_raw["Xs"][:,0:simulated_X.shape[0]]
full_or_est_X = V.numpy() @ simulated_X.T

err    = full_or_est_X - true_X
sq_err = (err**2).sum(0)
mse    = np.mean(sq_err)
print(f"\nMSE between simulated states and true states: {mse:.6f}")

y = np.zeros(len(test_DU)-1)
H = np.zeros(len(test_DU)-1)
for i in range(len(test_DU)-1):
    x    = torch.from_numpy(sol.y[:,i]).to(dtype=torch.float32)
    u    = test_DU[i].to(dtype=torch.float32)
    o1,y[i],H[i] = mor_PH.forward(x, u)

plt.rcParams.update({'font.size': 16})

plt.figure(figsize=(8, 5))
plt.plot(test_splt[:len(sq_err)], sq_err,  label="Squared Error of full-order state",    color="blue")
plt.xlabel('Time [s]')
# plt.yscale('log')
plt.legend()
plt.savefig(os.path.join("..","figs", "Error.pdf"))

plt.figure(figsize=(8, 5))
plt.plot(test_splt[:len(y)], test_DY[:len(y)], label = "Actual output",    color="red")
plt.plot(test_splt[:len(y)], y,                label = "Estimated output", color="blue", linestyle='--')
plt.xlabel('Time [s]')
plt.legend()
plt.savefig(os.path.join("..","figs", "Output.pdf"))

plt.figure(figsize=(8, 5))
plt.plot(test_splt[:len(H)], test_H[:len(H)], label = "Actual Hamiltonian", color="red")
plt.plot(test_splt[:len(H)], H, label = "Estimated Hamiltonian", color="blue", linestyle='--')
plt.xlabel('Time [s]')
plt.legend()
plt.savefig(os.path.join("..","figs", "Hamiltonian.pdf"))

r = V.shape[1]
def left(input):
    x   = torch.tensor(np.append(input,np.zeros(r-2)), dtype=torch.float32)
    x   = x.clone().detach().requires_grad_(True)
    NN  = mor_PH.H_net(x)
    dNN = autograd.grad(NN, x, create_graph=True)[0]
    return np.linalg.norm(dNN.detach().numpy(), ord=np.inf)

def right(input):
    x     = torch.tensor(np.append(input,np.zeros(r-2)), dtype=torch.float32)
    x_2   = math.sqrt(input[0]**2 + input[1]**2)
    x_inf = max(abs(input[0]), abs(input[1]))
    NN    = mor_PH.H_net(x).detach().numpy()
    return NN*d_smstep(input)*x_inf/(smstep(input)*math.sqrt(x_2**2 + sms_del**2))