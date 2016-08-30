import pkp.empirical_model
import matplotlib.pyplot as plt
import numpy as np

m = pkp.empirical_model.DAEM()
m.parameters = {'A0': 1e8, 'E0': 100e6, 'sigma': 5e6, 'y0': 0.6}

T = 1200
op_cond = [[0, T], [0.005, T]]

m.operating_conditions = op_cond


t, y = m.run()

plt.plot(t, y[:, 0])
# plt.show()

Em = m._Em
k = m.parameters['A0'] * np.exp(-Em / 8314.33 / T)
kint = 0
kint_v = np.empty((len(t), 4))
kint_v[0] = kint * k
for i, ti in enumerate(t[1:], 1):
    dt = ti - t[i - 1]
    kint += k * dt
    kint_v[i] = kint


coeff1 = m.Wm * m.mt / np.sqrt(np.pi)
coeff2 = np.exp(-np.power((Em - m.parameters['E0']) /
                          m.parameters['sigma'], 2) / 2)
coeff3 = np.exp(-kint_v)

ym = coeff1 * coeff2 * coeff3

y1 = 1 - np.sum(ym, axis=1)
y1 *= m.parameters['y0']

plt.plot(t, y1, marker='+', linewidth=0)
plt.show()
