import numpy as np
from pathlib import Path
from beamy import Support
from beamy.setup import plot_supports

# Create gallery directory
gallery_dir = Path("gallery")
gallery_dir.mkdir(exist_ok=True)

# Define beam length
L = 10.0

# Define various supports
supports = [
    Support(x=0.0, type="111111"),   # Fixed
    Support(x=3.0, type="111100"),   # Pinned
    Support(x=7.0, type="011000"),   # Roller
    # Support(x=L, type="000000")      # Free end (just a marker)
]

print("Generating support plot...")
plot_supports(
    supports, 
    beam_length=L,
    unit="m",
    save_path=str(gallery_dir / "supports_example.svg"),
    show=False
)

print(f"Plot saved to: {gallery_dir / 'supports_example.svg'}")

