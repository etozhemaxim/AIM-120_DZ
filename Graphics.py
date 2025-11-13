import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Чтение данных
data = pd.read_csv('all_angles_data.csv')

# Функция для расчета коэффициента торможения (аналогичная C++ версии)
def calculate_kappa_q_nos_con(M, lambda_nos):
    # Проверка граничных условий
    if lambda_nos < 0.7 or lambda_nos > 5.0:
        return 1.0
    
    # Случай дозвуковых скоростей
    if M <= 1.0:
        return 1.0
    
    # Данные для M=2.0
    lambda_values_M2 = [0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    kappa_values_M2 = [0.95, 0.92, 0.88, 0.85, 0.82, 0.80, 0.78, 0.76, 0.75, 0.74]
    
    # Интерполяция для M=2
    kappa_at_M2 = np.interp(lambda_nos, lambda_values_M2, kappa_values_M2)
    
    # Интерполяция между M=1 и M=2
    if M <= 2.0:
        return 1.0 + (M - 1.0) * (kappa_at_M2 - 1.0)
    else:
        # Для M > 2 используем значение при M=2
        return kappa_at_M2

# Параметры фюзеляжа (должны совпадать с C++ кодом)
l_nos = 0.47  # длина носовой части
D = 0.178     # диаметр
lambda_nos = l_nos / D

# Применяем торможение потока к данным
data_with_kappa = data.copy()

# Добавляем столбец с коэффициентом торможения
data_with_kappa['kappa_q'] = data_with_kappa['Mach'].apply(
    lambda M: calculate_kappa_q_nos_con(M, lambda_nos)
)

# Применяем торможение ко всем коэффициентам подъемной силы
angles = [-10, -5, 0, 5, 10]
for angle in angles:
    col = f'alpha_{angle}'
    data_with_kappa[col + '_with_kappa'] = data_with_kappa[col] * data_with_kappa['kappa_q']

# Создание фигуры с двумя графиками бок о бок
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# Цвета для графиков
colors = ['black', 'blue', 'green', 'orange', 'red', 'purple']
mach_colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'cyan', 'magenta']

# ЛЕВЫЙ ГРАФИК: C_y vs M для разных α (С ТОРМОЖЕНИЕМ И БЕЗ)
for i, angle in enumerate(angles):
    col = f'alpha_{angle}'
    
    # Без торможения (исходные данные) - пунктирная линия
    ax1.plot(data['Mach'], data[col], 
             color=colors[i % len(colors)], 
             linewidth=2,
             linestyle='--',
             alpha=0.7,
             label=fr'$\alpha = {angle}^\circ$ (без торможения)')
    
    # С торможением - сплошная линия
    ax1.plot(data_with_kappa['Mach'], data_with_kappa[col + '_with_kappa'], 
             color=colors[i % len(colors)], 
             linewidth=3,
             linestyle='-',
             label=fr'$\alpha = {angle}^\circ$ (с торможением)')

ax1.set_xlabel('M',fontsize=15 )
ax1.set_ylabel(r'$C_y\text{из.ф}$', fontsize=16)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper right', fontsize=10)

# ПРАВЫЙ ГРАФИК: C_y vs α для разных M (С ТОРМОЖЕНИЕМ И БЕЗ)
selected_mach = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

print(f"Выбранные числа Маха: {selected_mach}")

# Строим правый график
for j, mach in enumerate(selected_mach):
    # Находим ближайшее число Маха в данных
    idx = (data['Mach'] - mach).abs().idxmin()
    closest_mach = data['Mach'].iloc[idx]
    
    # Данные без торможения
    mach_data = data.iloc[idx]
    cy_values_no_kappa = []
    for angle in angles:
        col = f'alpha_{angle}'
        cy_values_no_kappa.append(mach_data[col])
    
    # Данные с торможением
    mach_data_with_kappa = data_with_kappa.iloc[idx]
    cy_values_with_kappa = []
    for angle in angles:
        col = f'alpha_{angle}_with_kappa'
        cy_values_with_kappa.append(mach_data_with_kappa[col])
    
    # Без торможения - пунктир
    ax2.plot(angles, cy_values_no_kappa, 
             color=mach_colors[j % len(mach_colors)], 
             linewidth=2, 
             linestyle='--',
             alpha=0.7,
             marker='s',
             markersize=4,
             label=fr'$M = {mach}$ (без торможения)')
    
    # С торможением - сплошная
    ax2.plot(angles, cy_values_with_kappa, 
             color=mach_colors[j % len(mach_colors)], 
             linewidth=3, 
             linestyle='-',
             marker='o',
             markersize=5,
             label=fr'$M = {mach}$ (с торможением)')

ax2.set_xlabel('$\\alpha$' ,fontsize=16)
ax2.set_ylabel(r'$C_y\text{из.ф}$', fontsize=16)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper left', fontsize=10) 

# Добавляем линию при α = 0 для наглядности
ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('comprehensive_aerodynamics_with_kappa.png', dpi=300, bbox_inches='tight')
plt.show()


# Чтение данных для изолированного крыла
data = pd.read_csv('krylo_isP.csv')

# Создание фигуры с двумя графиками бок о бок
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# Цвета для графиков
colors = ['black', 'blue', 'green', 'orange', 'red', 'purple']
mach_colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'cyan', 'magenta']

angles = [-10, -5, 0, 5, 10]

# ЛЕВЫЙ ГРАФИК: C_y vs M для разных α
for i, angle in enumerate(angles):
    col = f'alpha_{angle}'
    
    ax1.plot(data['Mach'], data[col], 
             color=colors[i % len(colors)], 
             linewidth=2,
             linestyle='-',
             label=fr'$\alpha = {angle}^\circ$')

ax1.set_xlabel('M', fontsize=15)
ax1.set_ylabel(r'$C_y\text{из.кр}$', fontsize=16)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper right', fontsize=10)

# ПРАВЫЙ ГРАФИК: C_y vs α для разных M
selected_mach = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

print(f"Выбранные числа Маха: {selected_mach}")

# Строим правый график
for j, mach in enumerate(selected_mach):
    # Находим ближайшее число Маха в данных
    idx = (data['Mach'] - mach).abs().idxmin()
    closest_mach = data['Mach'].iloc[idx]
    
    # Данные
    mach_data = data.iloc[idx]
    cy_values = []
    for angle in angles:
        col = f'alpha_{angle}'
        cy_values.append(mach_data[col])
    
    ax2.plot(angles, cy_values, 
             color=mach_colors[j % len(mach_colors)], 
             linewidth=2, 
             linestyle='-',
             marker='o',
             markersize=5,
             label=fr'$M = {mach}$')

ax2.set_xlabel('$\\alpha$', fontsize=16)
ax2.set_ylabel(r'$C_y\text{из.кр}$', fontsize=16)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper left', fontsize=10) 

# Добавляем линию при α = 0 для наглядности
ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('isolated_wing_aerodynamics.png', dpi=300, bbox_inches='tight')
plt.show()

# Чтение данных для изолированного крыла
data = pd.read_csv('krylo_isP_Intrf.csv')

# Создание фигуры с двумя графиками бок о бок
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# Цвета для графиков
colors = ['black', 'blue', 'green', 'orange', 'red', 'purple']
mach_colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'cyan', 'magenta']

angles = [-10, -5, 0, 5, 10]

# ЛЕВЫЙ ГРАФИК: C_y vs M для разных α
for i, angle in enumerate(angles):
    col = f'alpha_{angle}'
    
    ax1.plot(data['Mach'], data[col], 
             color=colors[i % len(colors)], 
             linewidth=2,
             linestyle='-',
             label=fr'$\alpha = {angle}^\circ$')

ax1.set_xlabel('M', fontsize=15)
ax1.set_ylabel(r'$C_y\text{из.кр}$', fontsize=16)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper right', fontsize=10)

# ПРАВЫЙ ГРАФИК: C_y vs α для разных M
selected_mach = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

print(f"Выбранные числа Маха: {selected_mach}")

# Строим правый график
for j, mach in enumerate(selected_mach):
    # Находим ближайшее число Маха в данных
    idx = (data['Mach'] - mach).abs().idxmin()
    closest_mach = data['Mach'].iloc[idx]
    
    # Данные
    mach_data = data.iloc[idx]
    cy_values = []
    for angle in angles:
        col = f'alpha_{angle}'
        cy_values.append(mach_data[col])
    
    ax2.plot(angles, cy_values, 
             color=mach_colors[j % len(mach_colors)], 
             linewidth=2, 
             linestyle='-',
             marker='o',
             markersize=5,
             label=fr'$M = {mach}$')

ax2.set_xlabel('$\\alpha$', fontsize=16)
ax2.set_ylabel(r'$C_y\text{из.кр}$', fontsize=16)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper left', fontsize=10) 

# Добавляем линию при α = 0 для наглядности
ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('isolated_wing_aerodynamics.png', dpi=300, bbox_inches='tight')
plt.show()



# Вывод информации о данных
print(f"Диапазон чисел Маха: от {data['Mach'].min()} до {data['Mach'].max()}")
print(f"Количество точек по M: {len(data)}")