import os
import numpy as np
import xarray as xr

# -------------------------------
# Utilities
# -------------------------------
def list_tmp_files(pthtmp):
    files = sorted([
        f for f in os.listdir(pthtmp)
        if f.startswith("RTOFS008_mesh025_gmapi_") and f.endswith(".npz")
    ], key=lambda x: int(x.split('_')[-1].split('.')[0]))
    return files


def load_tmp_files(pthtmp):
    IMOM, JMOM = [], []
    INDX_list, JNDX_list = [], []

    files = list_tmp_files(pthtmp)
    print(f"Found {len(files)} tmp files")

    for f in files:
        path = os.path.join(pthtmp, f)
        print(f"Loading {path}")
        data = np.load(path)

        IMOM.extend(data["imom"])
        JMOM.extend(data["jmom"])
        INDX_list.append(data["indx"])
        JNDX_list.append(data["jndx"])

    return IMOM, JMOM, INDX_list, JNDX_list


def save_tmp(dfltmp, IMOM, JMOM, INDX_list, JNDX_list):
    np.savez(
        dfltmp,
        imom=np.array(IMOM),
        jmom=np.array(JMOM),
        indx=np.concatenate(INDX_list, axis=0),
        jndx=np.concatenate(JNDX_list, axis=0)
    )


# -------------------------------
# Main processing loop
# -------------------------------
def process_domain(iS, iE, jS, jE,
                   hlon, hlat, HH,
                   LON, LAT,
                   regn, lat_min, lat_max, lat_sh, lat_nh,
                   use_tmp, dfltmp,
                   IMOM, JMOM, INDX_list, JNDX_list):

    pnts_saved = set(zip(IMOM, JMOM))
    n_failed = 0

    for ii in range(iS, iE):
        for jj in range(jS, jE + 1):

            if HH[jj, ii] >= 0:
                continue

            y0 = hlat[jj, ii]

            # region filtering
            if regn == 'global':
                if lat_nh > y0 > lat_sh:
                    continue
            else:
                if y0 < lat_min or y0 > lat_max:
                    continue

            if (ii, jj) in pnts_saved:
                continue

            x0 = hlon[jj, ii]

            ixx, jxx = mrmom.find_gridpnts_box(
                x0, y0, LON, LAT,
                dhstep=.8,
                ignore_north_lim=(lat_max >= 90.),
                wrap_long=True
            )
import os
import numpy as np
import xarray as xr

# -------------------------------
# Utilities
# -------------------------------
def list_tmp_files(pthtmp):
    files = sorted([
        f for f in os.listdir(pthtmp)
        if f.startswith("RTOFS008_mesh025_gmapi_") and f.endswith(".npz")
    ], key=lambda x: int(x.split('_')[-1].split('.')[0]))
    return files


def load_tmp_files(pthtmp):
    IMOM, JMOM = [], []
    INDX_list, JNDX_list = [], []

    files = list_tmp_files(pthtmp)
    print(f"Found {len(files)} tmp files")

    for f in files:
        path = os.path.join(pthtmp, f)
        print(f"Loading {path}")
        data = np.load(path)

        IMOM.extend(data["imom"])
        JMOM.extend(data["jmom"])
        INDX_list.append(data["indx"])
        JNDX_list.append(data["jndx"])

    return IMOM, JMOM, INDX_list, JNDX_list


def save_tmp(dfltmp, IMOM, JMOM, INDX_list, JNDX_list):
    np.savez(
        dfltmp,
        imom=np.array(IMOM),
        jmom=np.array(JMOM),
        indx=np.concatenate(INDX_list, axis=0),
        jndx=np.concatenate(JNDX_list, axis=0)
    )


# -------------------------------
# Main processing loop
# -------------------------------
def process_domain(iS, iE, jS, jE,
                   hlon, hlat, HH,
                   LON, LAT,
                   regn, lat_min, lat_max, lat_sh, lat_nh,
                   use_tmp, dfltmp,
                   IMOM, JMOM, INDX_list, JNDX_list):

    pnts_saved = set(zip(IMOM, JMOM))
    n_failed = 0

    for ii in range(iS, iE):
        for jj in range(jS, jE + 1):

            if HH[jj, ii] >= 0:
                continue

            y0 = hlat[jj, ii]

            # region filtering
            if regn == 'global':
                if lat_nh > y0 > lat_sh:
                    continue
            else:
                if y0 < lat_min or y0 > lat_max:
                    continue

            if (ii, jj) in pnts_saved:
                continue

            x0 = hlon[jj, ii]

            ixx, jxx = mrmom.find_gridpnts_box(
                x0, y0, LON, LAT,
                dhstep=.8,
                ignore_north_lim=(lat_max >= 90.),
                wrap_long=True
            )

            if len(ixx) == 0 or len(jxx) == 0:
                n_failed += 1
                continue

            ixx = np.expand_dims(ixx, axis=0)
            jxx = np.expand_dims(jxx, axis=0)

            IMOM.append(ii)
            JMOM.append(jj)
            INDX_list.append(ixx)
            JNDX_list.append(jxx)

            pnts_saved.add((ii, jj))   # ✅ important fix

    print(f"Discarded (no mapping): {n_failed}")

    if use_tmp:
        print(f"Saving tmp -> {dfltmp}")
        save_tmp(dfltmp, IMOM, JMOM, INDX_list, JNDX_list)

    return IMOM, JMOM, INDX_list, JNDX_list


# -------------------------------
# Final NetCDF writer
# -------------------------------
def save_netcdf(IMOM, JMOM, INDX_list, JNDX_list, LON, LAT, out_path):

    IMOM = np.array(IMOM)
    JMOM = np.array(JMOM)
    INDX = np.concatenate(INDX_list, axis=0)
    JNDX = np.concatenate(JNDX_list, axis=0)

    npnts = len(IMOM)
    jdim, idim = LON.shape

    ds = xr.Dataset({
        "mom_indx": (["npoints"], IMOM),
        "mom_jndx": (["npoints"], JMOM),
        "gmapi_i":  (["npoints", "nvert"], INDX),
        "gmapi_j":  (["npoints", "nvert"], JNDX),
        "longit":   (["jdim", "idim"], LON),
        "latit":    (["jdim", "idim"], LAT),
    })

    print(f"Saving NetCDF -> {out_path}")
    ds.to_netcdf(out_path, format='NETCDF4', engine='netcdf4')


# -------------------------------
# Driver logic
# -------------------------------
if final:
    # load all tmp files
    IMOM, JMOM, INDX_list, JNDX_list = load_tmp_files(pthtmp)

    # full domain pass
    iS, iE = 0, IDIM

    IMOM, JMOM, INDX_list, JNDX_list = process_domain(
        iS, iE, jS, jE,
        hlon, hlat, HH,
        LON, LAT,
        regn, lat_min, lat_max, lat_sh, lat_nh,
        use_tmp=False,
        dfltmp=None,
        IMOM=IMOM, JMOM=JMOM,
        INDX_list=INDX_list, JNDX_list=JNDX_list
    )

    save_netcdf(IMOM, JMOM, INDX_list, JNDX_list, LON, LAT, dfgmapi)

else:
    # chunk mode
    iS = (kchunk - 1) * chsize
    iE = kchunk * chsize

    IMOM, JMOM, INDX_list, JNDX_list = [], [], [], []

    if use_tmp and os.path.isfile(dfltmp):
        data = np.load(dfltmp)
        IMOM = data["imom"].tolist()
        JMOM = data["jmom"].tolist()
        INDX_list = [data["indx"]]
        JNDX_list = [data["jndx"]]

    process_domain(
        iS, iE, jS, jE,
        hlon, hlat, HH,
        LON, LAT,
        regn, lat_min, lat_max, lat_sh, lat_nh,
        use_tmp=True,
        dfltmp=dfltmp,
        IMOM=IMOM, JMOM=JMOM,
        INDX_list=INDX_list, JNDX_list=JNDX_list
    )
            if len(ixx) == 0 or len(jxx) == 0:
                n_failed += 1
                continue

            ixx = np.expand_dims(ixx, axis=0)
            jxx = np.expand_dims(jxx, axis=0)

            IMOM.append(ii)
            JMOM.append(jj)
            INDX_list.append(ixx)
            JNDX_list.append(jxx)

            pnts_saved.add((ii, jj))   # ✅ important fix

    print(f"Discarded (no mapping): {n_failed}")

    if use_tmp:
        print(f"Saving tmp -> {dfltmp}")
        save_tmp(dfltmp, IMOM, JMOM, INDX_list, JNDX_list)

    return IMOM, JMOM, INDX_list, JNDX_list


# -------------------------------
# Final NetCDF writer
# -------------------------------
def save_netcdf(IMOM, JMOM, INDX_list, JNDX_list, LON, LAT, out_path):

    IMOM = np.array(IMOM)
    JMOM = np.array(JMOM)
    INDX = np.concatenate(INDX_list, axis=0)
    JNDX = np.concatenate(JNDX_list, axis=0)

    npnts = len(IMOM)
    jdim, idim = LON.shape

    ds = xr.Dataset({
        "mom_indx": (["npoints"], IMOM),
        "mom_jndx": (["npoints"], JMOM),
        "gmapi_i":  (["npoints", "nvert"], INDX),
        "gmapi_j":  (["npoints", "nvert"], JNDX),
        "longit":   (["jdim", "idim"], LON),
        "latit":    (["jdim", "idim"], LAT),
    })

    print(f"Saving NetCDF -> {out_path}")
    ds.to_netcdf(out_path, format='NETCDF4', engine='netcdf4')


# -------------------------------
# Driver logic
# -------------------------------
if final:
    # load all tmp files
    IMOM, JMOM, INDX_list, JNDX_list = load_tmp_files(pthtmp)

    # full domain pass
    iS, iE = 0, IDIM

    IMOM, JMOM, INDX_list, JNDX_list = process_domain(
        iS, iE, jS, jE,
        hlon, hlat, HH,
        LON, LAT,
        regn, lat_min, lat_max, lat_sh, lat_nh,
        use_tmp=False,
        dfltmp=None,
        IMOM=IMOM, JMOM=JMOM,
        INDX_list=INDX_list, JNDX_list=JNDX_list
    )

    save_netcdf(IMOM, JMOM, INDX_list, JNDX_list, LON, LAT, dfgmapi)

else:
    # chunk mode
    iS = (kchunk - 1) * chsize
    iE = kchunk * chsize

    IMOM, JMOM, INDX_list, JNDX_list = [], [], [], []

    if use_tmp and os.path.isfile(dfltmp):
        data = np.load(dfltmp)
        IMOM = data["imom"].tolist()
        JMOM = data["jmom"].tolist()
        INDX_list = [data["indx"]]
        JNDX_list = [data["jndx"]]

    process_domain(
        iS, iE, jS, jE,
        hlon, hlat, HH,
        LON, LAT,
        regn, lat_min, lat_max, lat_sh, lat_nh,
        use_tmp=True,
        dfltmp=dfltmp,
        IMOM=IMOM, JMOM=JMOM,
        INDX_list=INDX_list, JNDX_list=JNDX_list
    )
