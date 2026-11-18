import numpy as np
import math
from scipy.integrate import quad


#Structure Variables
VerticalBeamLength = 52
VerticalBeamDiameter = 4

DiagonalBeamDiameter = 1

HorizontalBeamLength = 30 # from effective length
HorizontalBeamDiameter = 2

Depth = 40
d = Depth

Mass = 8*10**6 # kg

#Global Variables
g = 9.81 # m/s^2
Rho = 1025 # kg/m^3

#Waves Variables
WavePeriod = 10 # s
WaveHeight = 12 # m
WaveAmplitude = WaveHeight/2 # m
a = WaveAmplitude # m
MeanCurrent = 1.6 # m/s
omega = 2*np.pi/WavePeriod # rad/s



#Other Variables 
Cd = 1.2
Cm = 0.8

#-------------------------
#Part 1a - KC and maximum horizontal wave-induced drag and inertia forces

#Vertical Column
k_guess = omega**2 / g                  # deep-water initial guess
for i in range(100):
    k = omega**2 / (g * np.tanh(k_guess * d))
    if abs(k - k_guess) < 1e-10:    # converged
        break
    k_guess = k

print("Value of k: ", k)
print("Value of k * d: ", k * d)

waveLength = 2*np.pi/k
print("Value of wave length: ", waveLength)
print("--------------------------------")

#<CHECK!!> that we set z = 0 here
uMax = WaveAmplitude * omega * np.cosh(k * (d )) / np.sinh(k * d) * 1 # we set Z = 0 
print("Value of uMax: ", uMax) 

print("--------------------------------")
print("Part 1a i - Vertical Beam")

KCVerticalHorizontalMax = uMax * WavePeriod/ VerticalBeamDiameter 
print("Value of Vertical KC: ", KCVerticalHorizontalMax)

# Defining Eta, u and w
def eta(t):
    return(WaveAmplitude * np.cos(omega * t))

def u(z, t):
    return(WaveAmplitude * omega * np.cosh(k * (d + z)) / np.sinh(k * d) * np.cos(omega * t))

def w(z, t):
    return(WaveAmplitude * omega * np.sinh(k * (d + z)) / np.sinh(k * d) * np.sin(omega * t))

integral, _ = quad(lambda z: np.cosh(k * (z + d))**2, -d, 6.0)
Fd_max = 0.5 * Cd * Rho * VerticalBeamDiameter * WaveAmplitude**2 * omega**2 / np.sinh(k * d)**2 * integral

print("Value of Fd_max: ", Fd_max)


integral, _ = quad(lambda z: np.sinh(k * (z + d))**2, -d, 6.0)
Fi_max = Cm * Rho * ((np.pi*VerticalBeamDiameter**2 ) / (k*4)) *a*(omega**2)
print("Value of Fi_max: ", Fi_max)
#Values align with pics kelvin sent ;)

###########################################
print("--------------------------------\n")
print("Part 1a ii - Horizontal Beam")

HorizontalBeamDepth = -35

uMaxHorizontal = u(HorizontalBeamDepth, 0)

print("Value of u(HorizontalBeamDepth, 0): ", u(HorizontalBeamDepth, 0))

HorizontalKC = uMaxHorizontal * WavePeriod/ HorizontalBeamDiameter 

print("Value of Horizontal KC: ", HorizontalKC)

Fd_max = 0.5 * Cd * Rho * HorizontalBeamDiameter * HorizontalBeamLength * uMaxHorizontal**2
print("Value of Fd_max: ", Fd_max)



# Inertia: Cm rho (pi D^2 / 4) L * du/dt max at the beam depth
dudtMaxHorizontal = a * omega**2 * np.cosh(k * (d + HorizontalBeamDepth)) / np.sinh(k * d)
Fi_max = Cm * Rho * np.pi * HorizontalBeamDiameter**2 / 4 * HorizontalBeamLength * dudtMaxHorizontal

print("Value of Fi_max: ", Fi_max)
print("--------------------------------\n")

###########################################
print("Part 1a iii - Diagonal Beam ")

# Establishing geometric approximations
# Vertical separation = 47 - (1.5 + 1.5 + 5) = 39 m   (52 - 5 = 47 from the top)
# Horizontal separation = 30 - (2 + 2) = 26 m
VerticalSeparation = 47 - (1.5 + 1.5 + 5)
HorizontalSeparation = HorizontalBeamLength - VerticalBeamDiameter

DiagonalBeamLength = math.sqrt(HorizontalSeparation**2 + VerticalSeparation**2)
print("Value of DiagonalBeamLength: ", DiagonalBeamLength)

alpha = math.atan(VerticalSeparation / HorizontalSeparation)   # radians
sinAlpha = math.sin(alpha)
cosAlpha = math.cos(alpha)
print("Value of alpha (deg): ", math.degrees(alpha)) is 

# Diagonal runs from z = -33.5 (5 + 1.5 above the bed) to z = 5.5
zBottom = -d + 5 + 1.5
zTop = zBottom + VerticalSeparation

# KC = U T / D, where U is the velocity perpendicular to the beam axis
# Diagonals lie in the y-z plane, so u is fully perpendicular; w contributes w*cos(alpha)
uSWL = a * omega * np.cosh(k * d) / np.sinh(k * d)      # horizontal amplitude at z = 0
wSWL = a * omega * np.sinh(k * d) / np.sinh(k * d)      # vertical amplitude at z = 0 (= a*omega)
wPerpendicular = wSWL * cosAlpha
print("Value of u at SWL: ", uSWL)
print("Value of w at SWL: ", wSWL)
print("Value of perpendicular component of w: ", wPerpendicular)

UPerpMax = max(uSWL, wPerpendicular)
DiagonalKC = UPerpMax * WavePeriod / DiagonalBeamDiameter
print("Value of Diagonal KC: ", DiagonalKC)

# Drag: Fd = 1/2 Cd rho D a^2 w^2 / sinh^2(kd) * integral of cosh^2[k(z+d)] dz, z = -33.5 to 5.5
# Max at the wave crest (eta = 6 > 5.5), so the whole diagonal is submerged
# From Laplace tables: integral = [sinh(2k(z+d)) / 4k + z/2]
dragIntegral = (np.sinh(2 * k * (zTop + d)) - np.sinh(2 * k * (zBottom + d))) / (4 * k) + (zTop - zBottom) / 2
print("Value of drag integral: ", dragIntegral)

FdVerticalProjection = 0.5 * Cd * Rho * DiagonalBeamDiameter * a**2 * omega**2 / np.sinh(k * d)**2 * dragIntegral
print("Value of Fd over vertical projection: ", FdVerticalProjection)

# Integral was over dz; length along the beam is dz / sin(alpha)
Fd_max = FdVerticalProjection / sinAlpha
print("Value of Fd_max: ", Fd_max)

# Inertia: Fm = Cm rho (pi D^2 / 4) a w^2 / sinh(kd) * integral of cosh[k(z+d)] dz, z = -33.5 to 0
# du/dt is max at SWL (eta = 0), so the beam is only wet up to z = 0
inertiaIntegral = (np.sinh(k * (0 + d)) - np.sinh(k * (zBottom + d))) / k
print("Value of inertia integral: ", inertiaIntegral)

FmVerticalProjection = Cm * Rho * np.pi * DiagonalBeamDiameter**2 / 4 * a * omega**2 / np.sinh(k * d) * inertiaIntegral
print("Value of Fm over vertical projection: ", FmVerticalProjection)

Fm_max = FmVerticalProjection / sinAlpha
print("Value of Fm_max: ", Fm_max)

print("--------------------------------")

#-------------------------
#Part 1b - KC and maximum vertical wave-induced drag and inertia forces
# Vertical velocity:      w     = a w  sinh[k(z+d)] / sinh(kd) * sin(wt)
# Vertical acceleration:  dw/dt = a w^2 sinh[k(z+d)] / sinh(kd) * cos(wt)

###########################################
print("\nPart 1b i - Vertical Beam")

# w acts along the axis of the column. Morison drag and inertia only act
# perpendicular to a member's axis, so there is no vertical force on the column.
print("Vertical flow is parallel to the column axis: no vertical drag or inertia force")
print("--------------------------------\n")

###########################################
print("Part 1b ii - Horizontal Beam")

# w is fully perpendicular to the horizontal beam, evaluated at the beam depth z = -35
wMaxHorizontal = a * omega * np.sinh(k * (d + HorizontalBeamDepth)) / np.sinh(k * d)
print("Value of w max at beam depth: ", wMaxHorizontal)

HorizontalKCVertical = wMaxHorizontal * WavePeriod / HorizontalBeamDiameter
print("Value of Horizontal KC (vertical): ", HorizontalKCVertical)

# Drag: 1/2 Cd rho D L w^2
Fd_max_vertical = 0.5 * Cd * Rho * HorizontalBeamDiameter * HorizontalBeamLength * wMaxHorizontal**2
print("Value of vertical Fd_max: ", Fd_max_vertical)

# Inertia: Cm rho (pi D^2 / 4) L * dw/dt max at the beam depth
dwdtMaxHorizontal = a * omega**2 * np.sinh(k * (d + HorizontalBeamDepth)) / np.sinh(k * d)
Fi_max_vertical = Cm * Rho * np.pi * HorizontalBeamDiameter**2 / 4 * HorizontalBeamLength * dwdtMaxHorizontal
print("Value of vertical Fi_max: ", Fi_max_vertical)
print("--------------------------------\n")

###########################################
print("Part 1b iii - Diagonal Beam")

# Only the component of w perpendicular to the beam loads it: w*cos(alpha)
# KC uses that perpendicular component at SWL (found in part 1a iii)
DiagonalKCVertical = wPerpendicular * WavePeriod / DiagonalBeamDiameter
print("Value of perpendicular component of w at SWL: ", wPerpendicular)
print("Value of Diagonal KC (vertical): ", DiagonalKCVertical)

# Drag: force per unit length perpendicular to the beam = 1/2 Cd rho D (w cos(alpha))^2
# w is max when sin(wt) = 1, which is when eta = 0, so the beam is wet from z = -33.5 to 0
# integral of sinh^2[k(z+d)] dz = [sinh(2k(z+d)) / 4k - z/2]
dragIntegralVertical = (np.sinh(2 * k * (0 + d)) - np.sinh(2 * k * (zBottom + d))) / (4 * k) - (0 - zBottom) / 2
print("Value of vertical drag integral: ", dragIntegralVertical)

# Perpendicular force over the vertical projection, then / sin(alpha) for length along the beam
FdPerpendicular = 0.5 * Cd * Rho * DiagonalBeamDiameter * a**2 * omega**2 / np.sinh(k * d)**2 * cosAlpha**2 * dragIntegralVertical / sinAlpha
print("Value of Fd perpendicular to beam: ", FdPerpendicular)

# Vertical component of the perpendicular force
Fd_max_vertical = FdPerpendicular * cosAlpha
print("Value of vertical Fd_max: ", Fd_max_vertical)

# Inertia: force per unit length perpendicular to the beam = Cm rho (pi D^2 / 4) dw/dt cos(alpha)
# dw/dt is max when cos(wt) = 1, the wave crest (eta = 6 > 5.5), so the whole beam is submerged
# integral of sinh[k(z+d)] dz = [cosh(k(z+d)) / k], z = -33.5 to 5.5
inertiaIntegralVertical = (np.cosh(k * (zTop + d)) - np.cosh(k * (zBottom + d))) / k
print("Value of vertical inertia integral: ", inertiaIntegralVertical)

FmPerpendicular = Cm * Rho * np.pi * DiagonalBeamDiameter**2 / 4 * a * omega**2 / np.sinh(k * d) * cosAlpha * inertiaIntegralVertical / sinAlpha
print("Value of Fm perpendicular to beam: ", FmPerpendicular)

Fm_max_vertical = FmPerpendicular * cosAlpha
print("Value of vertical Fm_max: ", Fm_max_vertical)

print("--------------------------------")