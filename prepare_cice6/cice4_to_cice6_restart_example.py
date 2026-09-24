#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from netCDF4 import Dataset

# ==============================================================================
# CICE / Icepack Thermodynamic Constants
# ==============================================================================
L_SUB  = 2.835e6            # Latent heat of sublimation (J/kg)
L_VAP  = 2.501e6            # Latent heat of vaporization (J/kg)
LFRESH = L_SUB - L_VAP      # Latent heat of melting of fresh ice (J/kg)
RHOS   = 330.0              # Density of snow (kg/m^3)
CP_ICE = 2106.0             # Specific heat of fresh ice (J/kg/K)
RHO_W  = 1026.0             # Density of seawater (kg/m^3)
C_W    = 4218.0             # Specific heat of seawater (J/kg/K)
PUNY   = 1.0e-11            # Tiny volume threshold to prevent div-by-zero

def get_ice_salinity(z, smax=3.2, a=0.407, b=0.573):
    """Computes the Bitz and Lipscomb (1999) prescribed vertical salinity profile."""
    return 0.5 * smax * (1.0 - np.cos(np.pi * z**(a / (z + b))))

def clamp_mushy_enthalpy(qice, sice):
    """
    Caps the B&L enthalpy to prevent CICE6 Mushy-thermo warnings (zTin > Tmax).
    Calculates the mushy liquidus (Tmax) and caps enthalpy at qmlt - 1.0.
    """
    Tmax = np.zeros_like(sice)
    valid = sice > 0.0
    
    # Mushy liquidus curve
    Tmax[valid] = -sice[valid] / (18.48 - 0.01848 * sice[valid])
    
    # Enthalpy at melting point is purely sensible heat of the brine
    qmlt = RHO_W * C_W * Tmax
    
    # Cap the output enthalpy strictly below the melting point
    qice_clamped = np.copy(qice)
    qice_clamped[valid] = np.minimum(qice[valid], qmlt[valid] - 1.0)
    return qice_clamped

def read_record(fid, nx, ny):
    """Read a single 2D Fortran sequential binary record."""
    _ = np.fromfile(fid, dtype='>i4', count=1)
    data = np.fromfile(fid, dtype='>f8', count=nx*ny).reshape((ny, nx), order='C')
    _ = np.fromfile(fid, dtype='>i4', count=1)
    return data

def remap_enthalpy(q4, nilyr_in, nilyr_out):
    """Area-weighted vertical interpolation from input layers to output layers."""
    q6 = np.zeros((nilyr_out, *q4.shape[1:]), dtype=np.float64)
    for k in range(nilyr_out):
        top6, btm6 = k/nilyr_out, (k+1)/nilyr_out
        for j in range(nilyr_in):
            top4, btm4 = j/nilyr_in, (j+1)/nilyr_in
            overlap = max(0, min(btm6, btm4) - max(top6, top4))
            if overlap > 0:
                q6[k] += q4[j] * overlap * nilyr_out
    return q6

def main():
    parser = argparse.ArgumentParser(description="Convert RTOFS CICE4 binary restarts to CICE6 NetCDF format.")
    parser.add_argument("infile", help="Input CICE4 binary restart file")
    parser.add_argument("outfile", help="Output CICE6 NetCDF restart file")
    parser.add_argument("date", help="Target date in YYYYMMDD format")
    
    # Grid dimensions (defaults match 0.08-degree config)
    parser.add_argument("--nx", type=int, default=4500, help="Grid X dimension (default: 4500)")
    parser.add_argument("--ny", type=int, default=3297, help="Grid Y dimension (default: 3297)")
    parser.add_argument("--ncat", type=int, default=5, help="Number of ice categories (default: 5)")
    parser.add_argument("--nilyr_in", type=int, default=4, help="Number of ice layers in input (default: 4)")
    parser.add_argument("--nilyr_out", type=int, default=7, help="Number of ice layers in output (default: 7)")
    parser.add_argument("--nslyr", type=int, default=1, help="Number of snow layers (default: 1)")
    
    # Thermodynamics Option
    parser.add_argument("--ktherm", type=int, choices=[1, 2], default=1, 
                        help="Target thermodynamics: 1=Bitz&Lipscomb (default), 2=Mushy. Triggers enthalpy clamping if 2.")
    
    args = parser.parse_args()

    if len(args.date) != 8 or not args.date.isdigit():
        print("ERROR: Date must be exactly 8 digits (YYYYMMDD).")
        sys.exit(1)
        
    yr, mo, dy = int(args.date[0:4]), int(args.date[4:6]), int(args.date[6:8])
    spval = 1.e30

    print(f"Reading CICE4 binary restart: {args.infile}")
    with open(args.infile, 'rb') as fid:
        _ = np.fromfile(fid, dtype='>i4', count=1)
        istep = np.fromfile(fid, dtype='>i4', count=1)[0]
        _ = np.fromfile(fid, dtype='>f8', count=2)
        _ = np.fromfile(fid, dtype='>i4', count=1)

        aicen = np.array([read_record(fid, args.nx, args.ny) for _ in range(args.ncat)])
        vicen = np.array([read_record(fid, args.nx, args.ny) for _ in range(args.ncat)])
        vsnon = np.array([read_record(fid, args.nx, args.ny) for _ in range(args.ncat)])
        trcrn = np.array([read_record(fid, args.nx, args.ny) for _ in range(args.ncat)])

        eicen = np.array([read_record(fid, args.nx, args.ny) for _ in range(args.nilyr_in * args.ncat)])
        esnon = np.array([read_record(fid, args.nx, args.ny) for _ in range(args.nslyr * args.ncat)])

        uvel = read_record(fid, args.nx, args.ny)
        vvel = read_record(fid, args.nx, args.ny)
        scale_factor = read_record(fid, args.nx, args.ny)
        swvdr = read_record(fid, args.nx, args.ny)
        swvdf = read_record(fid, args.nx, args.ny)
        swidr = read_record(fid, args.nx, args.ny)
        swidf = read_record(fid, args.nx, args.ny)
        strocnxT = read_record(fid, args.nx, args.ny)
        strocnyT = read_record(fid, args.nx, args.ny)

        # CICE4 writes in order 1, 3, 2, 4
        stressp_1 = read_record(fid, args.nx, args.ny); stressp_3 = read_record(fid, args.nx, args.ny)
        stressp_2 = read_record(fid, args.nx, args.ny); stressp_4 = read_record(fid, args.nx, args.ny)
        
        stressm_1 = read_record(fid, args.nx, args.ny); stressm_3 = read_record(fid, args.nx, args.ny)
        stressm_2 = read_record(fid, args.nx, args.ny); stressm_4 = read_record(fid, args.nx, args.ny)
        
        stress12_1 = read_record(fid, args.nx, args.ny); stress12_3 = read_record(fid, args.nx, args.ny)
        stress12_2 = read_record(fid, args.nx, args.ny); stress12_4 = read_record(fid, args.nx, args.ny)

        iceumask = read_record(fid, args.nx, args.ny)
        _ = read_record(fid, args.nx, args.ny)
        _ = read_record(fid, args.nx, args.ny)

    print("Processing fields and mapping to CICE6...")
    
    def mask(arr): return np.where(arr > 0.5 * spval, 0.0, arr)
    
    uvel, vvel, scale_factor = mask(uvel), mask(vvel), mask(scale_factor)
    swvdr, swvdf, swidr, swidf = mask(swvdr), mask(swvdf), mask(swidr), mask(swidf)
    strocnxT, strocnyT, iceumask = mask(strocnxT), mask(strocnyT), mask(iceumask)
    
    stressp = [mask(stressp_1), mask(stressp_2), mask(stressp_3), mask(stressp_4)]
    stressm = [mask(stressm_1), mask(stressm_2), mask(stressm_3), mask(stressm_4)]
    stress12 = [mask(stress12_1), mask(stress12_2), mask(stress12_3), mask(stress12_4)]
    
    aicen, vicen, vsnon, trcrn = mask(aicen), mask(vicen), mask(vsnon), mask(trcrn)
    eicen, esnon = mask(eicen), mask(esnon)

    coszen = np.zeros((args.ny, args.nx), dtype=np.float64)
    fsnow = np.zeros((args.ny, args.nx), dtype=np.float64)
    iage = np.zeros_like(aicen)
    alvl = np.where(aicen > 0.0, 1.0, 0.0)
    vlvl = np.where(vicen > 0.0, 1.0, 0.0)
    apnd, hpnd, ipnd, dhs, ffrac = [np.zeros_like(aicen) for _ in range(5)]

    # Snow Enthalpy
    qsnon = np.zeros((args.nslyr, args.ncat, args.ny, args.nx))
    for n in range(args.ncat):
        vsn = np.maximum(vsnon[n], PUNY)
        qsn = (esnon[n] * args.nslyr) / vsn
        qsn = np.where(qsn > -RHOS * LFRESH, -RHOS * LFRESH, qsn)
        qsn = np.where(aicen[n] <= PUNY, 0.0, qsn)
        qsnon[0, n] = qsn

    # Ice Enthalpy (Input layers)
    qicen_in = np.zeros((args.nilyr_in, args.ncat, args.ny, args.nx))
    for n in range(args.ncat):
        vin = np.maximum(vicen[n], PUNY)
        for k in range(args.nilyr_in):
            qin = (eicen[n * args.nilyr_in + k] * args.nilyr_in) / vin
            qicen_in[k, n] = np.where(aicen[n] <= PUNY, 0.0, qin)

    # Remap Enthalpy 
    qicen_out = remap_enthalpy(qicen_in, args.nilyr_in, args.nilyr_out)

    # Ice Salinity Profile
    sice_out = np.zeros((args.nilyr_out, args.ncat, args.ny, args.nx))
    for k in range(1, args.nilyr_out + 1):
        z = (k - 0.5) / args.nilyr_out
        S = get_ice_salinity(z)
        for n in range(args.ncat):
            sice_out[k-1, n] = np.where(aicen[n] < 1.e-10, 0.0, S)
            
    # Conditional Thermo Safety Clamp
    if args.ktherm == 2:
        print("Applying Mushy thermodynamics enthalpy safety clamp...")
        qicen_out = clamp_mushy_enthalpy(qicen_out, sice_out)

    print(f"Writing CICE6 NetCDF restart: {args.outfile}")
    with Dataset(args.outfile, 'w', format='NETCDF4') as nc:
        nc.createDimension('ni', args.nx)
        nc.createDimension('nj', args.ny)
        nc.createDimension('ncat', args.ncat)

        nc.istep1 = istep
        nc.myear = yr
        nc.mmonth = mo
        nc.mday = dy
        nc.msec = 0

        def write_2d(name, data):
            var = nc.createVariable(name, 'f8', ('nj', 'ni'))
            var[:] = data

        def write_3d(name, data):
            var = nc.createVariable(name, 'f8', ('ncat', 'nj', 'ni'))
            var[:] = data

        write_2d('uvel', uvel); write_2d('vvel', vvel); write_2d('coszen', coszen)
        write_2d('scale_factor', scale_factor)
        write_2d('swvdr', swvdr); write_2d('swvdf', swvdf)
        write_2d('swidr', swidr); write_2d('swidf', swidf)
        write_2d('strocnxT', strocnxT); write_2d('strocnyT', strocnyT)
        
        for i in range(4): write_2d(f'stressp_{i+1}', stressp[i])
        for i in range(4): write_2d(f'stressm_{i+1}', stressm[i])
        for i in range(4): write_2d(f'stress12_{i+1}', stress12[i])

        write_2d('iceumask', iceumask); write_2d('fsnow', fsnow)

        write_3d('aicen', aicen); write_3d('vicen', vicen)
        write_3d('vsnon', vsnon); write_3d('Tsfcn', trcrn)
        write_3d('iage', iage); write_3d('alvl', alvl); write_3d('vlvl', vlvl)
        write_3d('apnd', apnd); write_3d('hpnd', hpnd); write_3d('ipnd', ipnd)
        write_3d('dhs', dhs); write_3d('ffrac', ffrac)

        for k in range(args.nilyr_out):
            write_3d(f'sice00{k+1}', sice_out[k])
            write_3d(f'qice00{k+1}', qicen_out[k])
        
        for k in range(args.nslyr):
            write_3d(f'qsno00{k+1}', qsnon[k])

    print("Success: CICE6 Restart generated.")

if __name__ == "__main__":
    main()
