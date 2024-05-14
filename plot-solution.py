import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Загрузка данных из CSV файла
dataframe = pd.read_csv('results.csv')
x = dataframe['x']
t = dataframe['t']
u = dataframe['u']

# Создание 3D графика
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(x, t, u, cmap='viridis')

# Настройка осей и меток
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u')
ax.set_title('График функции u(x,t)')

# Отображение графика
plt.show()
