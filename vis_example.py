import numpy as np
import matplotlib.pyplot as plt

def generate_banana_data():
    """Generates dummy data to simulate your experiment."""
    # 1. AP MICE (Fluid): Small random vectors (Noise)
    # 5 slices, ~50 random tiny vectors per slice
    ap_data = []
    for _ in range(5):
        noise = np.random.uniform(-0.1, 0.1, (50, 3)) 
        ap_data.append(noise)

    # 2. XCX MICE (Graft): Strong Ventral-Dorsal vectors
    # 5 slices. Slices 0 & 4 are weaker (edges), 1-3 are strong (center).
    aptr_data = []
    for i in range(5):
        # Create vectors pointing mostly UP (0, 1) and DOWN (0, -1) in 2D polar
        # We add some random jitter to make it look realistic
        num_fibers = 200
        base_y = np.random.choice([1, -1], num_fibers) # Up or Down
        jitter_x = np.random.normal(0, 0.1, num_fibers)
        jitter_y = np.random.normal(0, 0.1, num_fibers)
        
        vectors = np.zeros((num_fibers, 3))
        vectors[:, 0] = jitter_x       # X component (Left-Right)
        vectors[:, 2] = base_y + jitter_y # Z component (Ventral-Dorsal)
        
        # Scale amplitude: Middle slices are strong (0.8), edges weak (0.2)
        strength = 0.8 if i in [1, 2, 3] else 0.2
        vectors *= strength
        aptr_data.append(vectors)
        
    return ap_data, aptr_data

def plot_rose_summary(ap_data, aptr_data):
    fig, axes = plt.subplots(nrows=2, ncols=5, 
                             subplot_kw={'projection': 'polar'}, 
                             figsize=(15, 6))
    
    # Common radial limit to force comparison
    MAX_LIMIT = 50.0 

    # --- ROW 1: AP MICE (Fluid) ---
    for i, vecs in enumerate(ap_data):
        ax = axes[0, i]
        magnitudes = np.linalg.norm(vecs, axis=1)
        angles = np.arctan2(vecs[:, 2], vecs[:, 0]) # Z vs X
        
        # Weighted histogram
        bars, bin_edges = np.histogram(angles, bins=20, weights=magnitudes, range=(-np.pi, np.pi))
        
        # Plot
        width = (2*np.pi) / 20
        ax.bar(bin_edges[:-1], bars, width=width, color='gray', alpha=0.5)
        
        ax.set_ylim(0, MAX_LIMIT)
        ax.set_yticks([]) # Hide radial numbers
        ax.set_title(f"AP Slice {i+1}\n(Fluid/Noise)", fontsize=9)
        ax.set_theta_zero_location("N")

    # --- ROW 2: XCX MICE (Graft) ---
    for i, vecs in enumerate(aptr_data):
        ax = axes[1, i]
        magnitudes = np.linalg.norm(vecs, axis=1)
        angles = np.arctan2(vecs[:, 2], vecs[:, 0])
        
        bars, bin_edges = np.histogram(angles, bins=20, weights=magnitudes, range=(-np.pi, np.pi))
        
        width = (2*np.pi) / 20
        # Color the graft Blue to make it pop
        ax.bar(bin_edges[:-1], bars, width=width, color='blue', alpha=0.8)
        
        ax.set_ylim(0, MAX_LIMIT)
        ax.set_yticks([])
        ax.set_title(f"XCX Slice {i+1}\n(Graft)", fontsize=9, color='blue')
        ax.set_theta_zero_location("N")

    plt.suptitle("FOD Directionality: AP (Fluid) vs XCX (Graft)", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()

# Run it
ap, aptr = generate_banana_data()
plot_rose_summary(ap, aptr)