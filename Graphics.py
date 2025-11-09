import pandas as pd
import matplotlib.pyplot as plt

# Чтение данных из CSV
data = pd.read_csv('aerodynamics_data.csv')

# Построение графиков
plt.figure(figsize=(10, 6))
plt.plot(data['x'], data['y'], 'b-', linewidth=2, label='C_y^α')
plt.xlabel('Число Маха')
plt.ylabel('Производная коэффициента подъемной силы')
plt.title('Аэродинамические характеристики')
plt.grid(True)
plt.legend()
plt.savefig('aerodynamics_plot.png', dpi=300, bbox_inches='tight')
plt.show()