import pkp.biopolimi
import matplotlib.pyplot as plt

ua = pkp.biopolimi.bioS1.ultimate_analysis
print(ua)
pa = {'FC': 45.1,
      'VM': 50.6,
      'Ash': 4.3,
      'Moist': 19.0}

oc = [[0, 300], [0.01, 1000], [0.02, 1000]]

bio = pkp.biopolimi.BioPolimi(
    proximate_analysis=pa,
    ultimate_analysis=ua,
    pressure=101325,
    name='biomass')
bio.operating_conditions = oc

res = bio.run()

print(res.columns)
res = res.reset_index()

res.plot(x='Time, s', y='CELL')
plt.show()
