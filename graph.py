import matplotlib.pyplot as plt

# === DATA EXTRACTED FROM YOUR IMAGE ===
# Based on the visual height of the bars in the provided photo
years = [
    '2014', '2016', '2018', '2020', '2021',
    '2022', '2023', 'Early 2024', 'July 2024', '2025'
]
users = [35, 100, 200, 400, 500, 700, 800, 900, 950, 1000]

# Set figure size
plt.figure(figsize=(12, 7))

# Create the Bar Chart
# Using a blue color similar to the gradient in your image
bars = plt.bar(years, users, color='#3b82f6', width=0.6, zorder=3)

# Axis labels and Title
plt.ylabel('Number of monthly active users in millions', fontsize=12, color='gray')
plt.title('Telegram Number of Users, by Year', fontsize=16, fontweight='bold', pad=20)

# Add a horizontal grid (like in the photo) behind the bars
plt.grid(axis='y', linestyle='--', alpha=0.4, zorder=0)

# Clean up borders (remove top and right spines)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_color('lightgray')
plt.gca().spines['bottom'].set_color('lightgray')

# Add values on top of each bar for clarity
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2, height + 15,
             f'{height}', ha='center', va='bottom', fontsize=10, fontweight='bold', color='#333333')

# Rotate x-axis labels for better readability
plt.xticks(rotation=45, color='gray')
plt.yticks(color='gray')

# Show the plot
plt.tight_layout()
plt.show()