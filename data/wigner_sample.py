import sys
import numpy as np
from numpy import linalg as la
from scipy.spatial.transform import Rotation
import random

from utils import utils


def Laguerre(n, x):
    ## This function calculates laguerre polynomial
    ## L = n!/[(n-m)! * (m!)**2] = n*(n-1)*...*(n-m+1)/(m!)**2 = n/1**2 * (n-1)/2**2 *...*(n-m+1)/m**2, 0 <= m <= n

    L = 1 # L=1 when m=0   
    for m in range(1, n + 1):
        r = 1
        for mm in range(1, m + 1):
            r *= float(n - mm + 1) / mm**2
        L += (-1)**m * r * x**m

    return L


def Wignerfunc(mu, temp):
    ## This function generates random position Q and momenta P to find uptdate coifficents
    ## This function calls Laguerre to calculate the polynomial
    ## This function returns accepted Q and P

    #print('\nFreq: %s\n' % mu)
    max_pop = 0.9999
    ex = mu / (0.69503 * temp) #vibrational temperature: ex=h*c*mu/(kb*T), 0.69503 convert cm-1 to K
    pop= 0
    lvl_pop = []
    n = -1
    while True:
        n += 1
        pop += np.exp(-1 * ex * n) * (1 - np.exp(-1 * ex))
        lvl_pop.append(pop[0]) #Note here pop is a numpy array, thus pop[0] is the float number
        # Here is how I obtained this equation:
        # calculate partion function, fP=np.exp(ex*-0.5) /( 1 - np.exp(ex*-1) )
        # calculate population, pop=np.exp(-1*ex*(n+0.5))/fP
        #print('wignerfunction:%d %f %f %f %f'%(n,ex,np.exp(-1*ex*n)*(1-np.exp(-1*ex)),pop,max_pop))
        if pop >= max_pop:
            break
    while True:
        random_state=random.uniform(0,pop) # random generate a state
        n = -1
        for i in lvl_pop:                  # but population is not uniformly distributed over states
            n += 1
            if random_state <= i:          # find the lowest state that has more population than the random state
                break
        Q = random.uniform(0,1)*10.0-5.0
        P = random.uniform(0,1)*10.0-5.0
        rho2 = 2 * (Q**2 + P**2)
        W = (-1)**n * Laguerre(n, rho2) * np.exp(-0.5 * rho2)
        R=random.uniform(0, 1)
        #print('N: %d Q: %f P: %f W: %f R: %f' % (n,Q,P,W,R))
        if W > R and W < 1:
            #print('N: %d Q: %f P: %f Rho^2: %f W: %f R: %f' % (n,Q,P,rho2/2,W,R))

            break
    
    return float(Q), float(P)


def Wigner(sample):
    ## This function is based on SHARC wigner.py
    ## This function does Wigner sampling for structure and velocity
    ## This function calls Wiguerfunc to find update coefficient
    ## This function returns initial condition as [[atom x y z v(x) v(y) v(z)],...]

    temp = sample['temp']
    nfreq = sample['nfreq']
    freqs = sample['freqs']
    natom = sample['natom']
    atoms = sample['atoms']
    xyz = sample['xyz']
    vib = sample['vib']
    rmass = sample['rmass']
    amass = sample['amass']
    achrg = sample['achrg']

    mu_to_hartree    = 4.55633518e-6   # 1 cm-1  = h*c/Eh = 4.55633518e-6 au
    ma_to_amu        = 1822.88852      # 1 g/mol = 1/Na*me*1000 = 1822.88852 amu
    bohr_to_angstrom = 0.529177249     # 1 Bohr  = 0.529177249 Angstrom

    ## Constants:
    ##            meaning                     unit             conversion
    ## h          Planck constant             [kg*m^2/s]
    ## h_bar	  reduced Planck constant     [kg*m^2/s]       h_bar = h / (2*pi)
    ## a0         Bohr radii                  [m]              a0=sqrt(Eh/h_bar^2)
    ## Eh         Hartree energy              [kg*m^2/s^2]     Eh=h_bar^2/(me*a0^2)
    ## Na         Avogadro's number           
    ## ma         molar mass                  [g/mol]
    ## m          atomic mass                 [kg]
    ## me         electron mass               [kg]
    ## c          speed of light              [m/s]
    ## v          velocity                    [m/s]            v=Eh/me*vau
    ## vau        velocity in atomic unit     [Bohr/au]
    ## w          angluar frequency           [s-1]            w=2*pi*nu
    ## nu         frequency                   [s-1]            nu = mu * c
    ## mu         waveunmber                  [cm-1]           mu = nu / c
    ## vib        normal modes                [Bohr]           Molcas prints out non-weighted coordinates in Bohr
    ## n          mass-weighted normal modes  [Bohr*amu^0.5]   n=vib*sqrt(m)
    ## R/R0       cartesian coordinate        [Bohr]
    ## Q/P        dimensionless coordinates and momenta
    ##
    ## Math note:
    ##
    ## Convert wavenumber[cm-1]to energy[au] E=h*nu=h*c*mu=mu*(h*c/Eh)*Eh
    ## Convert molar mass[g/mol] to atomic mass unit m=ma/Na*me*1000
    ##
    ## Convert position q[m] to atomic distance[Bohr]
    ## Q = sqrt(w*m/h_bar)*q => q = Q*sqrt(h_bar/w*m) = Q*sqrt(h_bar^2/w*h_bar*m) = Q*sqrt(h_bar^2/Eh*me)/[ sqrt(E/Eh)*sqrt(m/me) ] = Q*a0/[ sqrt(mu*(h*c/Eh))*sqrt(ma/(Na*me*1000)) ]
    ## q = Q/[ sqrt(mu*(h*c/Eh))*sqrt(ma/(Na*me*1000)) ] in [Bohr]
    ## Update position R from R0
    ## R = R0+sum(sum(q*n)*mode)
    ## q*n = Q/[ sqrt(mu*(h*c/Eh))*sqrt(ma/(Na*me*1000)) ]*vib*sqrt(m) = Q/sqrt(mu*(h*c/Eh))*sqrt(Na*me*1000)*vib
    ## q*n = Q/sqrt(mu*mu_to_hartree)*sqrt(1/ma_to_amu)*vib
    ##
    ## Convert velocity[m/s] to atomic unit[Bohr/au]
    ## P = p/sqrt(w*m*h_bar) => p = m*v = P*sqrt(w*m*h_bar) => v = P*sqrt(w*h_bar/m) = P*sqrt(Eh/me)*sqrt(E/Eh)/sqrt(m/me) = P*vau*sqrt(mu*(h*c/Eh))/sqrt(ma/Na*me*1000))
    ## v = P*sqrt(mu*(h*c/Eh))/sqrt(ma/Na*me*1000)) in [Bohr/au]
    ## Update velcocity from 0 to v 
    ## v = 0+sum(sum(p*n)*mode)
    ## v*n = P*sqrt*(mu*(h*c/Eh))/sqrt(ma/(Na*me*1000))*vib*sqrt(m) = P*sqrt(mu*(h*c/Eh))*sqrt(Na*me*1000)*vib
    ## v*n = P*sqrt(mu*mu_to_hartree)*sqrt(1/ma_to_amu)*vib

    Q_P = np.array([Wignerfunc(i, temp) for i in freqs])   # generates update coordinates and momenta pairs Q and P

    Q = Q_P[:, 0].reshape((nfreq, 1))                      # first column is Q

    Q *= 1 / np.sqrt(freqs * mu_to_hartree * ma_to_amu)    # convert coordinates from m to Bohr
    Qvib = np.array([np.ones((natom, 3)) * i for i in Q])  # generate identity array to expand Q
    Qvib = np.sum(vib * Qvib, axis = 0) * 0.3                    # sum sampled structure over all modes
    newc = (xyz + Qvib) * bohr_to_angstrom                 # cartesian coordinates in Angstrom

    # P = Q_P[:, 1].reshape((nfreq, 1))                      # second column is P
    # P *= np.sqrt(freqs * mu_to_hartree / ma_to_amu)        # convert velocity from m/s to Bohr/au
    # Pvib = np.array([np.ones((natom,3)) * i for i in P])   # generate identity array to expand P
    # velo = np.sum(vib * Pvib, axis = 0)                    # sum sampled velocity over all modes in Bohr/au
    #                                                        # use velo=1 in &DYNAMIX to read un-weighted velocity
    # #These lines are for test
    # #amass=np.array([np.ones(3)*i for i in amass])         # expand atoic mass over x y z coordinates
    # #print(np.linalg.norm(vib[0]*amass**0.5))              # check if mode is mass-weighted (norm=1)  
    # #Epot=0.5*np.sum((freqs*mu_to_hartree*Q)**2*ma_to_amu) # The Q has been updated, so it has to be converted back
    # #Ekin=0.5*np.sum(amass*ma_to_amu*velo**2)         
    # #print(Epot,Ekin,Epot+Ekin)

    # #velo*=np.sqrt(amass*ma_to_amu)                        # only use for velo=2 in &DYNAMIX convert velocity[Bohr/au] to mass-weighted [Bohr*amu^0.5/au] for Molcas

    # inicond = np.concatenate((newc,velo), axis = 1)
    # inicond = np.concatenate((atoms.reshape((-1, 1)), inicond), axis = 1)
    # inicond = np.concatenate((inicond, amass), axis = 1)
    # inicond = np.concatenate((inicond, achrg), axis = 1)

    return newc



#####################################
#
# MAP THE ATOMS
# F-C-C-F
# 0 1 2 3
# 
# Hs all at the end
# H1, H2 correspond to C1
# H3, H4 correspond to C2

C1_IDX = 1
C2_IDX = 2
F1_IDX = 0
F2_IDX = 3
H1_IDX = 4
H2_IDX = 5
H3_IDX = 6
H4_IDX = 7

#####################################

def rotate_difluoroethane(geom, n):

    bond_axis = geom[C1_IDX] - geom[C2_IDX]
    bond_axis_norm = np.linalg.norm(bond_axis)
    n_bond_axis = bond_axis / bond_axis_norm


    to_rotate = [F2_IDX, H3_IDX, H4_IDX]

    raw_angles = np.linspace(0, 360, n, endpoint=False)
    angles = []

    for idx, a in enumerate(raw_angles):

        if idx == 0:
            angles.append(a)

        elif idx == n-1:
            angles.append(a)

        else:
            previous = raw_angles[idx-1]
            next = raw_angles[idx+1]

            angles.append(np.random.uniform(previous, next))


    rotated = []

    for angle in angles:
        rotation = Rotation.from_rotvec(np.radians(angle) * n_bond_axis)

        updated = geom.copy()
        for idx in to_rotate:
            vec = geom[idx] - geom[C1_IDX]
            updated[idx] = rotation.apply(vec) + geom[C1_IDX]

        rotated.append(updated)

    return rotated


def stretch_difluoroethane(geom, n):

    bond_length = np.linalg.norm(geom[C1_IDX] - geom[C2_IDX])

    disp = bond_length * 0.05

    to_move = [C2_IDX, F2_IDX, H3_IDX, H4_IDX]

    raw_bond_lengths = np.linspace(bond_length - disp, bond_length + disp, n)
    bond_lengths = []

    for idx, l in enumerate(raw_bond_lengths):

        if idx == 0:
            bond_lengths.append(l)

        elif idx == n-1:
            bond_lengths.append(l)

        else:
            previous = raw_bond_lengths[idx-1]
            next = raw_bond_lengths[idx+1]

            bond_lengths.append(np.random.uniform(previous, next))

    stretched = []

    for bl in bond_lengths:

        scale = bl / bond_length
        updated = geom.copy()

        for idx in to_move:
            vec = geom[idx] - geom[C1_IDX]
            updated[idx] = geom[C1_IDX] + vec * scale

        stretched.append(updated)

    return stretched


def Sample(stem, n):

    info = utils.read_mopac_freq_info(stem)
    info["temp"] = 298.15

    ensemble = []

    for s in range(n):

        try:
            ensemble.append(Wigner(info))

        except Exception as e:
            print(str(e))

    return ensemble



if __name__ == "__main__":

    xyz = sys.argv[1]

    elements, coords = utils.read_xyz(xyz)
    stem = ".".join(xyz.split(".")[:-1])
    n = 200

    wigner_structures = Sample(stem, n)

    all = []

    for geom in wigner_structures:
        all.append(geom)

        for rotated in rotate_difluoroethane(geom, 50):
            all.append(rotated)

            for stretched in stretch_difluoroethane(rotated, 50):
                all.append(stretched)

    to_keep = int(min(10000, len(all)))
    to_keep = 10
    all = [all[i] for i in np.random.choice(list(range(len(all))), to_keep)]

    for idx, geom in enumerate(all):
        utils.write_xyz(elements, geom, f"sampled_structrue_{idx}.xyz")


