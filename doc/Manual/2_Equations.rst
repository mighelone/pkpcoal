Equations
=========

Pyrolysis Kinetic Equations
---------------------------

In PKP are four possible kinetic models to select:

-  The Constant Rate Model, section [SSS\_cR]

-  The Arrhenius Model, [SSS\_Arrh]

-  The Kobayashi Model, [SSS\_Kob]

-  The Distributed Activation Energy Model, [SSS\_DAEM]

| The Kobayashi model has the advantage, that the yields are dependent
  on the individual temperature history. A change, may in in the final
  temperature leads to a higher yield.
| The constant rate and the Arrhenius model has fixed final yields. The
  general devolatilization reactions for these two models are described
  with equations [E\_Reaction\_s] and [E\_Reaction\_g],
  where [E\_Reaction\_s] describes the mass loss of the
  solid (index \ *s*) coal, equation [E\_Reaction\_g] the formation of
  the yields (gaseous:\ *g*, individual:\ *i*).

  .. math::

     \begin{aligned}
     \label{E_Reaction_s}
      \frac{dm_{s}}{dt}&=-k_s \: \left( m_{s} - m_{s,final} \right)\\
     \label{E_Reaction_g}
      \frac{dm_{g,i}}{dt}&=k_{g,i} \: \left(m_{g,i,final} - m_{g,i}\right)  \end{aligned}

  So when using more than one run, also the final yield is a parameter
  to fit. It will be located in between the final yields of the runs of
  the pyrolysis models.

In the following subsections, the plots showing the rates and yields
were generated with and showing the yields and rates for different
species. All are based on the same temperature profiles as shown in
figure [F\_Tt].

.. figure:: Figures/tempHist
   :alt: The temperature history.

   The temperature history.

The table [T\_Fit] gives an overview of the models with their parameters
to fit and species to apply this model on. This is described more
detailed in the next subsections.

[T\_Fit]

+-----------------+-----------------------+----------------------+
| **Model**       | **fitted parmeter**   | **fitted species**   |
+=================+=======================+======================+
| constant Rate   |                       |                      |
+-----------------+-----------------------+----------------------+
| Arrhenius       |                       |                      |
+-----------------+-----------------------+----------------------+
| Kobyashi        |                       |                      |
+-----------------+-----------------------+----------------------+
| DAEM            |                       |                      |
+-----------------+-----------------------+----------------------+

Table: The fitted models and their parameter.

The Constant Rate Model
~~~~~~~~~~~~~~~~~~~~~~~

Assuming a **constant rate** (:math:`\mathrm{k = const. }` and a
starting time :math:`\mathrm{t_{start}}`), the
equations [E\_Reaction\_s] and [E\_Reaction\_g] can be solved
analytically:

.. math::

   \begin{aligned}
   \label{E_constRate_s}
   m_s(t)&=m_{s,final} + \left( m_{s}(t=t_{start,s}) - m_{s,final} \right) \: e^{-k_s(t-t_{start,s})}\\
   \label{E_constRate_g}
   m_{g,i}(t)&=m_{g,i,final}\cdot \left( 1 - \: e^{-k_i(t-t_{start,i})} \right)\\
   \label{E_Offset_Time}
   if \;\;\; t\leq t_i\::\;\;\; m(t)&=m(0)\end{aligned}

| This leads to advantage concerning the computational costs. On the
  other hand is this model completely independent from the temperature
  history. This is visualized in figure [F\_Fit\_cR\_Y] where all fitted
  curves overlap each other.
| The parameters to fit are for every species:

-  the starting time (\ :math:`t_{start,i}`), where the devolatilization
   begins

-  the constant rate factor :math:`k_i`

-  the final yield (when more than one run)

This model is applied to all species.

.. figure:: Figures/FG-DVC-Fit_result_cR_H2O_Y
   :alt: One fitting result (Yields) for the constant Rate Model. The
   observed species is water.

   One fitting result (Yields) for the constant Rate Model. The observed
   species is water.

The Arrhenius Model
~~~~~~~~~~~~~~~~~~~

| The kinetic rate k can be also expressed with the **Arrhenius**
  equation:

  .. math::

     \begin{aligned}
     \label{E_Arrhenius_s}
      \frac{dm_s}{dt}&=A_s \cdot T(t)^{b_s} \cdot e^{-\frac{E_s}{T(t)}}\left( m_{s} - m_{s,final} \right)\\
     \label{E_Arrhenius_g}
      \frac{dm_{g,i}}{dt}&=A_i \cdot T(t)^{b_{g,i}} \cdot e^{-\frac{E_{g,i}}{T(t)}}\left(m_{g,i,final} - m_{g,i}\right)\end{aligned}

  This notation of the Arrhenius equation includes no gas constant R in
  the exponential term. So the activation energy (or here more precise
  activation temperature) has the unit Kelvin. This is so far an
  advantage as the fitted :math:`E_i` is independent of the used unit
  system (SI or cgs).
| A second notation of the Arrhenius equation contains not the
  correction term :math:`T(t)^{b_{g,i}}`:

  .. math::

     \begin{aligned}
     \label{E_Arrhenius_s_noB}
      \frac{dm_s}{dt}&=A_s \cdot e^{-\frac{E_s}{T(t)}}\left( m_{s} - m_{s,final} \right)\\
     \label{E_Arrhenius_g_noB}
      \frac{dm_{g,i}}{dt}&=A_i \cdot e^{-\frac{E_{g,i}}{T(t)}}\left(m_{g,i,final} - m_{g,i}\right)\end{aligned}

Unlike the constant rate model is the Arrhenius modeled rate influenced
by the temperature. But the Arrhenius model can be used to express the
evolve for all species and the final yields are also fixed. So the
parameter to fit are here:

-  the preexponentiation factor \ :math:`A_i`

-  the correction factor :math:`b_i`

-  the activation energy :math:`E_i`

-  the final yield (when more than one run)

| This model is applied to all species.
| The Arrhenius model leads to a good agreement in the yield and rate
  curves for a limited range of temperatures,
  figures [F\_Fit\_Arrh\_Y], [F\_Fit\_Arrh\_R]. The disadvantage is the
  temperature independent yield fraction, all integrals for the rate
  curves in figure [F\_Fit\_Arrh\_R] are the same. This leads to an
  imprecision as the yields show a dependency on the final temperature
  and heating rate.

.. figure:: Figures/CPD-Fit_result_Arrh_Tar_Y
   :alt: One fitting result (Yields) for the Arrhenius Model. The
   observed species is tar.

   One fitting result (Yields) for the Arrhenius Model. The observed
   species is tar.

.. figure:: Figures/CPD-Fit_result_Arrh_Tar_R
   :alt: One fitting result (Rates) for the Arrhenius Model. The
   observed species is tar.

   One fitting result (Rates) for the Arrhenius Model. The observed
   species is tar.

The Kobayashi Model
~~~~~~~~~~~~~~~~~~~

Also the **Kobayashi** equation, also Two Competing Reaction Model, can
be fitted, see equation [E\_Kobayashi]. The optimization is carried out
using the Arrhenius notation of equation [E\_Kob\_k] for
:math:`\mathrm{k_1}` and :math:`\mathrm{k_2}`.

.. math::

   \label{E_Kobayashi}
    \frac{m_v(t)}{m_{p,0} - m_a}= \int_{0}^{t} ( \alpha_1 k_1 + \alpha_2 k_2 ) exp \left( -  \int_{0}^{t} ( k_1 + k_2 ) \; dt \right) \; dt

.. math::

   \label{E_Kob_k}
    k_j=A_j \:e^{-\frac{E_{j}}{T(t)}} \;\;\;\;\;\; with \: j=1,2

| The Kobayashi model can be applied only on the the overall, the total
  yields. The yield of individual species could be generated by
  multiplying the overall yield with the yield
  fraction \ :math:`\frac{y_i}{y_{all}}`. But as the composition of the
  yields varies with the temperature history this factor also shows this
  dependency, which may leads to an imprecision when modeling the
  individual yields.
| The final yields of this model are dependent on the temperature, see
  figures [F\_Fit\_Kob\_Y] and [F\_Fit\_Kob\_R]. The range of the yields
  are defined by the two weight
  factors \ :math:`\alpha_1` and :math:`\alpha_2`. The :math:`k_1`
  models the reactions at lower temperatures (low :math:`A_1` and
  :math:`E_1`), :math:`k_2` at higher temperatures (high :math:`A_2` and
  :math:`E_2`). If the final temperature has very low values, the yields
  will converge to \ :math:`\alpha_1`. If the temperatures will raise to
  infinity, the yields will be equal \ :math:`\alpha_2`. So
  :math:`\alpha_2` is ever set equal one: :math:`\alpha_2=1`.
  :math:`\alpha_1` is equal the amount of volatile matter in the daf
  coal. As the measurements, the proximate analysis of coal is based on,
  were carried out at very low heating rates compared with the ones
  occurring at gasification and combustion processes, the approximation
  :math:`\alpha_1=VM` is an applicable and good assumption.
| For the fitting procedure, the inner integral
  :math:`\int_{0}^{t} ( k_1 + k_2 ) \; dt` is approximated by the
  Trapezoidal rule.
| As it can be seen in figures [F\_Fit\_Kob\_Y] and [F\_Fit\_Kob\_R],
  the higher temperatures (figure [F\_Tt]) lead to higher yields. But
  the influence of the temperature cannot be modeled that the dependency
  on temperature is completely the same as in the output of the more
  complex pyrolysis programs. So leads the higher temperature in case 1
  compared to case 3 to a slightly higher yield in the output of while
  the influence on the Kobayashi modeled result is greater.

.. figure:: Figures/CPD-Fit_result_Kob_Total_Y
   :alt: One fitting result (Yields) for the Kobayashi Model. The
   Kobayashi model just optimizes the overall yields.

   One fitting result (Yields) for the Kobayashi Model. The Kobayashi
   model just optimizes the overall yields.

.. figure:: Figures/CPD-Fit_result_Kob_Total_R
   :alt: One fitting result (Rates) for the Kobayashi Model.

   One fitting result (Rates) for the Kobayashi Model.

The Distributed Activation Energy Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| The Distributed Activation Energy Model (DAEM) considers parallel
  first order kinetics over a specific range, described by a
  distribution function (F(E)). The equation used in PKP is:

  .. math::

     \label{E_DAEM}
      m = m_{final} \left( 1 - \int_{0}^{\infty} exp\left[ -A_0 \cdot \int^{t_{final}}_{t_0} exp\left( -\frac{E}{T} \right) dt  \right] F(E) \right)

  As a distribution function, a Gaussian Distribution is used:

  .. math::

     \label{E_GaussDistr}
      F(E) = \frac{1}{\sigma \cdot \sqrt{2\pi}} \cdot exp \left( -\frac{(E-E_0)^2}{2\sigma^2} \right)

   So there are four parameters to optimize:

#. the preexponentiation factor :math:`\mathrm{A_0}`, which is here
   equal for all reactions.

#. :math:`\mathrm{E_0}`, defining the center of the Gaussian
   Distribution

#. and :math:`\mathrm{\sigma}`, which spcifies the flattening of the
   Distribution curve and its range

#. for multiple runs, the :math:`\mathrm{m_{final}}` also has to be
   optimized

| For solving the outer integral over dE, the Simpson tule is used. But
  for the nemuerical solution of the integral, the range has to be
  modified. As reported in the paper by
  Cai :raw-latex:`\cite{Cai_DAEM1}`, the integration boundaries can be
  se to :math:`\mathrm{E_0 + 3 \cdot \sigma}` as the upper and
  :math:`\mathrm{E_0 - 3 \cdot \sigma}` as the lower limit. This range
  covers up to 99.73 % of the applied Gaussian Distribution. So for the
  further fitting of the devolatilization reaction, the prerequisite
  :math:`\mathrm{E_0 > 3 \cdot \sigma}` should be in force to achieve
  realistic results.
| The inner integral
  :math:`\mathrm{\int^{t_{final}}_{t_0} exp\left( -\frac{E}{T} \right) dt}`
  was already simplyfied by setting :math:`\mathrm{A_0}` as a constant.
  In many
  papers :raw-latex:`\cite{Cai_DAEM1,Cai_DAEM2,Cai_DAEM3,Slovak_DAEM}` a
  linear heating rate over the whole time range is assumed, so
  :math:`\mathrm{\frac{dT}{dt}=\beta}`. The transformation of the
  integral allows to integrate over the temperature. For such
  temperature integrals different analytical approaches
  exist :raw-latex:`\cite{Cai_DAEM1,Cai_DAEM2}`. This method was also
  tested. But as this is only a very specific case of the operating
  conditions and even not faster than solving the integrals [1]_, this
  approach is furthernot considered. The double integral is solved
  numerically. [2]_
| The outer integral is solved over a specific number of activation
  energies. For each activation energy the inner integral is solved and
  the value for all time steps saved, using a compisite Trapezoidal
  rule [3]_. All inner integrals are saved in a
  2D-Array(\ :math:`\mathrm{t_i,E_i}`). Each column contains all values
  of the inner integral for all time steps (the same activation energy),
  each line allinner integrals at the time :math:`\mathrm{t_i}`. After
  this array is calculated, the equation of the outer integral is used
  for all values of the current line of the matrix and afterwards this
  list is integrated over dE.

Fitting the Kinetic Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fitting procedure is carried out with a ``scipy-optimizer``\  [4]_
 [5]_ and the ``scipy.odeint``\  [6]_ package to minimize the
residual \ :math:`E(k,m_{fit})` in the equation [E\_LS]. For the
structure of the whole fitting procedure see chapter [S\_Program]. In
equation [E\_LS] is :math:`\mathrm{m_{out}}` the output of the
devolatilization program or . The optimization is carried out over all
points reported in the output file of the pyrolysis program. The
normalized weight factor
parameters \ :math:`\mathrm{a_0}` and :math:`\mathrm{a_1}` in the
equations [E\_Weight\_Param1] and [E\_Weight\_Param2] can both be
defined by the user, the standard setting is for both one.

.. math::

   \label{E_LS}
    E(k,m_{fit})=\omega_0 \int \left( m_{out}(t) - m_{fit}(k,t) \right)^2 dt \; + \; \omega_1 \int \left( \dot m_{out}(t) - \dot m_{fit}(k,t) \right)^2 dt

.. math::

   \begin{aligned}
    \label{E_Weight_Param1}
   \omega_0 &= \frac{a_0}{\left( max(m_{out})-min(m_{out}) \right)^2}\\
    \label{E_Weight_Param2}
   \omega_1 &= \frac{a_1}{max(\dot{m}_{out}^2)}\end{aligned}

Pyrolysis Species- and Energy Conservation for output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Species Conservation
^^^^^^^^^^^^^^^^^^^^

As the first step it has to be checked that the oxygen content in the
generated yields (oxygen containing species :math:`\mathrm{f_i}` with
:math:`\mu_i^{O}` oxygen) is less equal that the oxygen in the ultimate
analysis:

.. math::

   UAO^{cpd, \: species \: output} = M_{O} \sum_i \frac{\mu_i^{O} f_i}{M_i} \le UAO 
    \label{E_O_balance}

The factor \ :math:`\mathrm{\gamma}` ([E\_gamma]) tells if the outputted
yields contain less oxygen than reported in the Ultimate
Analysis (UA) (\ :math:`\gamma > 1`) or if they are
equal (\ :math:`\gamma = 1`). For the case that \ :math:`\gamma < 1`,
the oxygen containing yields have to be decreased by using
equation [E\_scale\_up], while the amount of the other species have to
be increased to conserve the conserve the amount of volatile
matter (equation [E\_add\_up]). In this case, the tar will contain no
oxygen. The yield of :math:`N_2` is equal to the UA of Nitrogen.

.. math::

   \begin{aligned}
    \gamma &= \frac{UAO}{UAO^{cpd, \: species \: output}}
    \label{E_gamma}\\
    f_i^{new} &= \gamma \cdot f_i 
    \label{E_scale_up}\\
    f_{oth}^{new} &= f_{oth} + \left(1-\gamma\right) \sum_i f_i
    \label{E_add_up}\end{aligned}

For the case :math:`f_{N_2}<f_{other}`, the remaining part is assigned
to :math:`CH_4`:

.. math::

   f_{CH_4}^{new} = f_{CH_4} + \left( f_{oth}^{new} - f_{N_2} \right)
   \label{E_MethanNew}

Now the composition of tar can be calculated. For each element C,H,O,
the following equation [E\_TarComp] can be used, assuming a tar
composition of \ :math:`C_nH_mO_p`. :math:`M_j` is the atom weight of
the element j, :math:`\mu_i^j` the number of atoms of \ :math:`j` in the
species \ :math:`i`.

.. math::

   \frac{UA_j}{M_j} = \mu_{tar}^j \frac{f_{tar}}{M_{tar}} + \sum_i \mu_i^{j} \frac{f_{i}}{M_{i}}
   \label{E_TarComp}

Energy Conservation
^^^^^^^^^^^^^^^^^^^

| The Dulong formula is used, if the higher heating value (HHV) of the
  coal is not known:

  .. math::

     \label{E_Dulong}
      HHV = 32.79 \cdot UAC + 150.4 \cdot (UAH - UAO/8) + 9.26 \cdot UAS + 4.97 \cdot UAO + 2.42 \cdot UAN

  where UAC, UAH, UAO, UAS and UAN are the value of the ultimate
  analysis for carbon, hydrogen, oxygen, sulfur and nitrogen. The result
  has the unit of \ :math:`\frac{MJ}{kg_{coal, as recieved}}`.
| Afterwards, the HHV (entered by the user or calculated with the Dulong
  formula) for the coal as received is related to the dry ash-free (daf)
  state (equation [E\_HHVdaf]). This new HHV is used to get the lower
  heating value for a daf state, equation [E\_LHV]. In this equation,
  :math:`r_{H_2O}` is the latent heat of water.

  .. math::

     \begin{aligned}
      HHV_{daf}&=\frac{HHV_{ar}}{PAVM+PAFC}
     \label{E_HHVdaf}\\
     LHV_{daf}&=HHV_{daf}-\frac{M_{H_2O}}{2 \cdot M_H} \cdot UAH \cdot r_{H_2O}\
     \label{E_LHV}\end{aligned}

| Regarding the combustion of the raw coal (equation [E\_Raw\_Comb]),
  the energy balance can be written as in equation [E\_Raw\_hf].

  .. math::

     \begin{aligned}
      &C_xH_y O_z N_w + (x + y/4 - z/2) O2 \rightarrow x CO2 + y/2 H2O + w/2 N2
     \label{E_Raw_Comb}\\
     &Q_{react}=LHV_{raw}\cdot M_{daf} = h_{f,raw} + (x + y/4 - z/2) h_{f,O_2} -x h_{f,CO_2} -y/2
         h_{f,H_2O} - w/2 h_{f,N_2}
     \label{E_Raw_hf}\end{aligned}

  Using equation [E\_Raw\_hf], the heat of formation of the raw
  molecule (\ :math:`h_{f,raw}`) can be calculated.
| The heat of formation for tar is based on the equation [E\_DevolTar],
  implying, that no heat is produced or absorbed during the
  devolatilization process.

  .. math::

     \label{E_DevolTar}
      C_x H_y O_z N_w \rightarrow \nu_{char}C_{(s)} + \nu_{tar} C_n H_m O_p + \sum_i \nu_i M_i

| The stoichiometric coefficient of each species can be calculated from
  the volatile yield expressed as mass fraction:

  .. math::

     \label{E_myTar}
      \nu_i = \frac{f_i M_{raw}}{M_i}

  Making the energy balance for equation [E\_DevolTar] with
  :math:`Q_{react}=0`, the heat of formation for tar is:

  .. math::

     \label{E_Tar_hf}
      \nu_{tar} h_{f,tar} = h_{f,raw} - \nu_{char} h_{f,char} - \sum_i \nu_i h_{f,i}
| Another method is to assume a heat of formation for tar equal zero
  (e.g. if there is only a very low yield of tar), and calculate the
  heat of pyrolysis:

  .. math:: - Q_{pyro} \cdot M_{raw} = h_{f,raw} - \nu_{char} h_{f,char} - \sum_i \nu_i h_{f,i}

  Where :math:`Q_{pyro}` is the heat of pyrolysis per unit of mass of
  daf. It is positive if heat is required for breaking coal structure
  bounds. Generally, it is expressed in terms of volatile matter:

  .. math::

     \label{E_QPyro}
      Q_{pyro}^{vm} = \frac{Q_{pyro}}{1-f_{char}}

Pyrolysis Species and Energy Conservation for output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Species Conservation
^^^^^^^^^^^^^^^^^^^^

| As in most of the CFD applications the combustion of some output
  species like HCN, COS or Olefins are not implemented. Only the species
  Char, Tar, CO, :math:`CO_2`, :math:`H_2O`, :math:`CH_4` and
  :math:`H_2` are further considered. The nitrogen is merged into the
  tar. So the amount of tar is calculated by the using
  equation [E\_newTar], where the sum contains the species Char, CO,
  :math:`CO_2`, :math:`H_2O`, :math:`CH_4` and :math:`H_2`.

  .. math::

     \label{E_newTar}
      f_{Tar}=1-\sum_i f_i

  Applying the equation [E\_TarComp] for all the elements j (Carbon,
  Hydrogen, Nitrogen, Oxygen), the composition of tar (its
  stoichiometric coefficients) is calculated.

Energy Conservation
^^^^^^^^^^^^^^^^^^^

Applying the energy balance on the combustion reaction of the
devolatilization yields (for the case of a non-heat producing/consuming
pyrolysis process), the following reaction equation is satisfied:

.. math::

   \label{E_TarEnergy}
    LHV_{daf}=H_{f,Tar} \cdot f_{Tar} + \sum_i H_{f,i} \cdot f_{i}

The LHV is calculated based on the HHV using the same equations as in
chapter [SSS\_ConsEqCPD]. The :math:`H_{f,i}` are calculated with the
following equations [E\_hf1] to [E\_hf4], making an energy balance for
every of the pyrolysis yields. So the heat of formation for tar can be
calculated from equation [E\_TarEnergy], as all other parameters are
known.

.. math::

   \begin{aligned}
   \label{E_hf1}
    H_{f,Char}&=\left( (h_{f,Char}+h_{f,O_2}-h_{f,CO_2}) \cdot f_{Char} \right) \cdot M_C^{-1} \\
   \label{E_hf2}
    H_{f,H_2}&=\left( (h_{f,H_2}+ \frac{1}{2} \cdot h_{f,O_2} - h_{f,H_2O}) \cdot f_{H_2} \right) \cdot M_{H_2O}^{-1} \\
   \label{E_hf3}
    H_{f,CH_4}&=\left( (h_{f,CH_4}+ 2 \cdot h_{f,O_2}-h_{f,CO_2}-2 \cdot h_{f,H_2O}) \cdot f_{CH_4} \right) \cdot M_{CH_4}^{-1} \\
   \label{E_hf4}
    H_{f,CO}&=\left( (h_{f,CO}+ \frac{1}{2} \cdot h_{f,O_2}-h_{f,CO_2}) \cdot f_{CO} \right) \cdot M_{CO}^{-1}\end{aligned}

| The :math:`H_{f,Tar}` with the unit \ :math:`\frac{J}{kg}` is
  transformed back into \ :math:`\frac{J}{kmol}` by multiplying with the
  molecular mass of tar.
| To calculate the heat of formation for tar, the tar combustion can be
  regarded, as the tar composition is known:

  .. math:: C_nH_mO_pN_k + \nu_{O_2} O_2 \rightarrow  \nu_{CO_2} CO_2 + \nu_{H_2O} H_2O + \nu_{N_2} N_2

This leads to the balance, where the :math:`h_{f,Tar}` can be
calculated:

.. math:: H_{f,Tar} = h_{f,Tar} + \nu_{O_2} h_{f,O_2} - \nu_{CO_2} h_{f,CO_2} - \nu_{H_2O} h_{f,H_2O} -\nu_{N_2} h_{f,N_2}

.. [1]
   The reason might be that the very large analytical equations in
   Python (as an interpreting language) take more time than let the
   integral solve by a external library.

.. [2]
   There are also other approaches to avoid the double integration like
   in the paper by McGuiness et al. :raw-latex:`\cite{McGuiness_DAEM}`,
   but also these assumptions made here
   (:math:`\mathrm{\sigma \rightarrow 0 }` or
   :math:`\mathrm{\frac{E_0}{T} \rightarrow \infty }`) cannot be applied
   here.

.. [3]
   the ``scipy.integrate.cumtrapz`` module,
   http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.cumtrapz.html

.. [4]
   *http://docs.scipy.org/doc/scipy/reference/optimize.html*

.. [5]
   The standard setting optimizer is ``fmin``. The ``leastsq`` optimizer
   is the second choice.

.. [6]
   *http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html*
