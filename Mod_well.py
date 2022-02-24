#Etape01

def rayon_adimensionnel(re, rw):
  """Calculer le rayon adimensionnel (rD)"""
  return re / rw

def temps_adimensionnel(perm, t, poro, mu, ct, rw):
  """calculer le temps adimensionnel (tD)"""
  return (.0002637 * perm * t) / (poro * mu * ct * (rw**2))

def time_finite_acting(perm, poro, mu, ct, rw, re):
  """Calculer le moment où le flux commence à se comporter de manière infinie"""
  r_D = re / rw
  t_Dw = 0.25 * r_D**2
  return (poro * mu * ct * (rw**2) * t_Dw) / (.0002637 * perm)

def pression_adimensionnelle(rD, tD):
    """
    Calculer la pression sans dimension à partir d'un débit constant
    """
    import numpy as np
    if tD < (0.25 * rD**2):
        # Solution pour le cas d'action infinie pour un débit constant (Eq. 4.9)
        pD = 0.5 * (np.log(tD) + .80907)
    if tD > (0.25 * rD**2):
        # Solution pour le cas d'action finie pour un débit constant (Eq. 4.8)
        pD = (2 * tD / rD**2) + np.log(rD) - .75
    return pD

def debit_adimensionnel(rD, tD):
    """
    Calculer un débit adimensionnel à partir d'une pression constante
    """
    import numpy as np
    import pandas as pd
    from scipy.interpolate import griddata
    
    if tD < (0.25 * rD**2):
        # Solution pour le cas d'action infinie pour un débit constant (Eq. 4.22, 4.23)
        if tD > 0.01 and tD < 200:
          # Eq. 6.42
          qD = (26.7544 + (45.5537 * np.sqrt(tD)) + (13.3813 * tD) + (0.492949 * tD * np.sqrt(tD))) / ((47.4210 * np.sqrt(tD)) + (35.5372 * tD) + (2.60967 * tD * np.sqrt(tD)))
        if tD >= 200:
          # Eq. 6.43
          qD = ((2.02623 * tD * (np.log(tD) - 1)) + 3.90086) / (tD * ((np.log(tD))**2))

    if tD > (0.25 * rD**2):
        # Solution pour le cas d'action finie pour un débit constant
        qD = np.nan
        qD = 2 / (np.log(tD) + .80907)

    return qD
#---------------------------------------
#Etape02

def pression_multirate(pD, delta_q, pi, B, mu, perm, h):
  """Calculer la pression d'écoulement comme étant la somme de débits constants"""
  import numpy as np
  return pi - ((B * mu / (.007082 * perm * h)) * (np.sum(pD * delta_q)))

def debit_multipressure(qD, delta_p, B, mu, perm, h):
  """Calculer le débit comme étant la somme des pressions d'écoulement constantes"""
  import numpy as np
  return ((.007082 * perm * h) / (B * mu)) * (np.sum(qD * delta_p))
#---------------------------------------
#Etape03
def simulation_multirate_test(p_initial, t_step, t_change, q_change,
                            re, rw, perm, poro, mu, ct, Bo, h):
  """
  Simulation de l'essai à taux constant multiple 
  Basé sur le principe de superposition
  """
  import numpy as np
  import matplotlib.pyplot as plt
  import matplotlib.patches as mpl_patches
  
  # calcul du temps d'action fini
  t_finite_acting = time_finite_acting(perm, poro, mu, ct, rw, re)

  # produire le vecteur de temps
  t_end = t_change[-1]
  time = np.arange(0, t_end+1, t_step)

  # calculer le rayon adimensionnel
  rD = re / rw

  # calculer la différence de débit (Δq)
  t_change = np.append(0, t_change)
  delta_q = [j-i for i, j in zip(q_change[:-1], q_change[1:])]
  delta_q = np.concatenate((np.array([0, q_change[0]]), delta_q))

  # créer un profil de pas de taux
  tmax = t_change[-1] + 1
  t = []
  q = []
  pwf = []

  for i in range(len(time)):  
      for j in range(0, len(t_change)-1):
          if time[i] > t_change[j] and time[i] <= t_change[j+1]:
              # produire le profil t et q
              t.append(time[i])
              q.append(q_change[j])
              
              # calculer le temps adimensionnel tD (tD1, tD2, ..., tDn) à chaque fois
              tn = time[i] - t_change[:j+1]    
              tD = temps_adimensionnel(perm, tn, poro, mu, ct, rw)
              
              # calculer la pression adimensionnelle pD à chaque fois
              pD = []
              for k in range(len(tD)):
                  _ = pression_adimensionnelle(rD, tD[k])
                  # _ = pd(rD, tD[k])
                  pD.append(_)
              
              # calculer la pression finale après superposition
              delta_qn = delta_q[1:j+2] 
              
              pwf_ = pression_multirate(pD, delta_qn, p_initial, Bo, mu, perm, h)
              pwf.append(pwf_)       

  t, q, pwf = np.append(0, t), np.append(q_change[0], q), np.append(p_initial, pwf)

  # tracer le débit du puits et le profil de pression d'écoulement
  plt.figure(figsize=(17,5))

  ## sortir le temps d'action fini dans le graphe
  labels = []
  labels.append("Time @ Finite-acting = {} heurs".format(np.round(t_finite_acting, 2)))

  handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                  lw=0, alpha=0)] * 1

  ## tracer le débit
  plt.subplot(1,2,1)
  plt.step(t, q, color='blue')
  plt.title('Profil du débit de puits', size=20, pad=15)
  plt.xlim(0, t_end)
  plt.ylim(ymax=max(q)+200)
  plt.xlabel('Temps (heurs)'); plt.ylabel('Débit (STB/D)')

  plt.legend(handles, labels, loc='upper right', fontsize=12, 
              fancybox=True, framealpha=0.7, 
              handlelength=0, handletextpad=0) 

  ## tracer BHFP
  plt.subplot(1,2,2)
  # t = np.arange(len(pwf))
  plt.plot(t, pwf, color='red')
  plt.title('Profil de pression découlement du puits', size=20, pad=15)
  plt.xlim(0, t_end)
  plt.xlabel('Temps (heurs)'); plt.ylabel('BHFP (psia)')

  plt.show()
  
def simulate_multipressure_test(p_initial, t_step, t_change, p_change,
                                re, rw, perm, poro, mu, ct, Bo, h):
  """
  Simulation de l'essai de pression d'écoulement constante multiple (BHFP) 
  Basé sur le principe de superposition
  """
  import numpy as np
  import matplotlib.pyplot as plt
  import matplotlib.patches as mpl_patches
  
  # calcul du temps d'action fini
  t_finite_acting = time_finite_acting(perm, poro, mu, ct, rw, re)  

  # produire le vecteur de temps
  t_end = t_change[-1]
  time = np.arange(0, t_end+1, t_step)

  # calculer le rayon adimensionnel
  rD = re / rw

  # calculer la différence de débit (Δq)
  t_change = np.append(0, t_change)
  pi_min_p0 = p_initial - p_change[0]
  delta_p = [i-j for i, j in zip(p_change[:-1], p_change[1:])]
  delta_p = np.concatenate((np.array([0, pi_min_p0]), delta_p))

  # créer un profil de pas de taux
  tmax = t_change[-1] + 1
  t = []
  pwf = []
  q = []

  for i in range(len(time)):  
      for j in range(0, len(t_change)-1):
          if time[i] > t_change[j] and time[i] <= t_change[j+1]:
              # produire le profil t et p 
              t.append(time[i])
              pwf.append(p_change[j])
              
              # calculer le temps adimensionnel tD (tD1, tD2, ..., tDn) à chaque fois
              tn = time[i] - t_change[:j+1] # is an array   
              tD = temps_adimensionnel(perm, tn, poro, mu, ct, rw)
              
              # calculate le débit adimensionnel qD à chaque instant
              qD = []
              for k in range(len(tD)):
                  _ = debit_adimensionnel(rD, tD[k])
                  # _ = qd(rD, tD[k])
                  qD.append(_)
              
              # calculer le débit final après superposition
              delta_pn = delta_p[1:j+2] 
              
              q_ = debit_multipressure(qD, delta_pn, Bo, mu, perm, h)
              q.append(q_)     

  # tracer le profil de la pression d'écoulement et du débit du puits
  plt.figure(figsize=(17,5))

  ## sortir le temps d'action fini dans le graphe
  labels = []
  labels.append("Time @ Finite-acting = {} hours".format(np.round(t_finite_acting, 2)))

  handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                  lw=0, alpha=0)] * 1

  ## tracer BHFP
  plt.subplot(1,2,1)
  plt.step(t, pwf, color='red')
  plt.title('Profil de pression découlement du puits', size=20, pad=15)
  plt.xlim(0, t_end)
  plt.ylim(ymax=max(pwf)+200)
  plt.xlabel('Temps (heurs)'); plt.ylabel('Pression (psia)')

  plt.legend(handles, labels, loc='upper right', fontsize=12, 
              fancybox=True, framealpha=0.7, 
              handlelength=0, handletextpad=0) 

  ## tracer le débit
  plt.subplot(1,2,2)
  plt.plot(t, q, color='blue')
  plt.title('Profil du débit de puits', size=20, pad=15)
  plt.xlim(0, t_end)
  plt.xlabel('Temps (heurs)'); plt.ylabel('Débit (STB/D)')

  plt.show()