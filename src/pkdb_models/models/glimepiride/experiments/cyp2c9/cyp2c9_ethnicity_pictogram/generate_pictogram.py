import os
from PIL import Image
from pkdb_models.models.glimepiride import RESULTS_PATH
from pathlib import Path

icon_path = Path(__file__).resolve().parent / "person.png"

# Frequency data for each population (*1/*1, *1/*2, *1/*3, *3/*3, Other)
population_frequencies = {
    "Afr. Amer./Afro-Carib.": [75.87, 3.91, 2.35, 0.02, 17.85],
    "American": [83.15, 6.09, 5.48, 0.09, 5.19],
    "Central/South Asian": [59.62, 17.57, 16.96, 1.21, 4.64],
    "East Asian": [83.79, 0.39, 6.89, 0.14, 8.80],
    "European": [62.85, 20.18, 11.98, 0.57, 4.42],
    "Latino": [74.34, 13.15, 6.92, 0.16, 5.42],
    "Near Eastern": [61.13, 20.29, 12.90, 0.68, 5.00],
    "Oceanian": [91.22, 5.60, 2.98, 0.02, 0.18],
    "Sub-Saharan African": [52.64, 1.90, 1.62, 0.01, 43.83],
}

# Color for each genotype
cyp2c9_colors = {
    "*1/*1": "black",
    "*1/*2": "#b589d6",
    "*1/*3": "#804fb3",
    "*3/*3": "#552586",
    "Other": "grey",
}

output_dir = RESULTS_PATH / "cyp2c9_pictograms"
os.makedirs(output_dir, exist_ok=True)

# Grid setup for arranging icons
grid_rows, grid_cols = 5, 20
icon_size = 500

# Load and resize the person icon
icon = Image.open(icon_path).convert("RGBA").resize((icon_size, icon_size))

# Generate pictograms per population
for population, frequencies in population_frequencies.items():
    # Compute rounded count of each genotype
    counts = {genotype: round(frequencies[k]) for k, genotype in enumerate(cyp2c9_colors.keys())}

    # Flatten into a list with multiple entries per genotype
    icon_genotypes = []
    for genotype, count in counts.items():
        icon_genotypes += [genotype] * count

    # Setup canvas size
    row_spacing, col_spacing = 5, 3
    canvas_width = grid_cols * (icon_size + col_spacing)
    canvas_height = grid_rows * (icon_size + row_spacing)
    canvas = Image.new("RGBA", (canvas_width, canvas_height), (255, 255, 255, 0))

    # Paste colored icons to canvas
    for index, genotype in enumerate(icon_genotypes):
        row = index // grid_cols
        col = index % grid_cols
        x = col * (icon_size + col_spacing)
        y = row * (icon_size + row_spacing)

        # Tint the icon according to genotype color
        tint = Image.new("RGBA", icon.size, cyp2c9_colors[genotype])
        colored_icon = Image.composite(tint, icon, icon)

        canvas.paste(colored_icon, (col * icon_size, row * icon_size), mask=colored_icon)

    # Resize and save each individual pictogram
    final_size = (400, 300)
    canvas = canvas.resize(final_size, resample=Image.LANCZOS)
    safe_name = population.replace("/", "_").replace(" ", "_").replace(".", "")
    canvas.save(os.path.join(output_dir, f"{safe_name}.png"))

# Combine all individual pictograms into one composite image
pictogram_files = sorted([
    f for f in output_dir.glob("*.png")
    if not f.stem.startswith("cyp2c9_ethnicities_combined")
])
images = [Image.open(f).convert("RGBA") for f in pictogram_files]

# Use first image size as base for layout
img_width, img_height = images[0].size
spacing = 20

# Calculate final width of combined image
total_width = len(images) * img_width + (len(images) - 1) * spacing
composite = Image.new("RGBA", (total_width, img_height), (255, 255, 255, 255))

# Paste all individual images side by side
for index, img in enumerate(images):
    x = index * (img_width + spacing)
    composite.paste(img, (x, 0), mask=img)

# Save the final combined image
composite.save(output_dir / "cyp2c9_ethnicities_combined.png")