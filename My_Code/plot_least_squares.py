import pandas as pd
import matplotlib.pyplot as plt
import os

# Get the absolute path to the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
# Construct the full path to the CSV file
csv_path = os.path.join(script_dir, 'least_squares_data.csv')

data = pd.read_csv(csv_path)
plt.figure(figsize=(10, 6))
plt.plot(data['x'], data['y'], 'k*', markersize=10, label='Data points')
plt.plot(data['x'], data['predicted_y'], 'r-', linewidth=2, label='LS1 (optimal)')
plt.plot(data['x'], data['alt_predicted_y'], 'g-.', linewidth=2, label='LS2 (linalg)')
plt.plot(data['x'], data['not_as_good_y'], 'b--', linewidth=2, label='Suboptimal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Simple Least Squares: y = b*x')
plt.grid(True)
plt.legend()
plt.savefig('least_squares_plot.png')
plt.show()
