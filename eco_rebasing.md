# Rebasing ECO code to `lalsuite-v6.73`

*(c) Yasmeen Asali, July 2020*  

## Setup

First step, setting up a new `lalsuite` installation on CIT. Get an SSH key for gitlab:
```bash
ssh-keygen -o -t rsa -b 4096 -C "yasmeen.asali@LIGO.ORG"
pbcopy < ~/.ssh/id_rsa.pub
```
Paste contents [here](https://git.ligo.org/profile/keys). 

Create a virtual environment: 
```bash
mkdir -p Projects/eco_resonance
cd Projects/eco_resonance
virtualenv -p python3.6 ecoENV/
unset PYTHONPATH
source ecoENV/bin/activate
pip install pip setuptools -U
pip install numpy scipy matplotlib healpy h5py ipython
```

Clone `lalsuite`:
```bash 
git clone git@git.ligo.org:eco-resonance-tgr/lalsuite.git
cd lalsuite/
git checkout -b eco-resonance origin/eco-resonance
```

Check that I'm on the correct branch by listing all branches:
```bash
git branch -a
```

Make `lalsuite`:
```bash
./00boot
./configure --prefix=/home/yasmeen.asali/Projects/eco_resonance/ecoENV/opt/lalsuite --enable-python --enable-swig-python
make -j
make install
```
Use `-j` to make with all available cores. 

The first time you `make install` you will need to add the following to the end of the `ecoENV/bin/activate` script:
```bash 
. /home/yasmeen.asali/Projects/eco_resonance/ecoENV/opt/lalsuite/etc/lalsuite-user-env.sh
```
This workd for bourne shells. For other shells, follow the instructions printed at the end of `make install`.

Check that `lalsuite` is working properly in ipython:
```python
import lal
import lalsimulation as lalsim

m1 = 25*lal.MSUN_SI
m2 = 20*lal.MSUN_SI

#check rerquired arguments
lalsim.SimIMRPhenomDGenerateFD?

#generate waveform
output = lalsim.SimIMRPhenomDGenerateFD(0., 20., 0.25, m1, 
                                        m2, 0., 0., 20., 1024., 100e6*lal.PC_SI)

htilde = output.data.data
```


## Adding ECO changes to LALSimulation

Four files to change: 
- `lalsimulation/lib/LALSimInspiral.c`
- `lalsimulation/lib/LALSimInspiral.h`
- `lalsimulation/lib/LALSimInspiralWaveformParams.c`
- `lalsimulation/lib/LALSimInspiralWaveformParams.h`

See Anurdha's changes [here](https://git.ligo.org/eco-resonance-tgr/lalsuite/-/commit/98d658027d0b850fd6bacc271d19b5d271e1dd61)

1 place in LALSimInspiral.c:
```c
// around line ~2273
    if (XLALSimInspiralWaveformParamsLookupEnableECO(LALparams))
      ret = XLALSimECOResonanceTerm(hptilde, hctilde, m1/LAL_MSUN_SI, m2/LAL_MSUN_SI, distance, LALparams);

// around line ~6232
int XLALSimECOResonanceTerm(
                                          COMPLEX16FrequencySeries **hptilde, /**< Frequency-domain waveform h+ */
                                          COMPLEX16FrequencySeries **hctilde, /**< Frequency-domain waveform hx */
                                          REAL8 m1,                           /**< Mass 1 in solar masses */
                                          REAL8 m2,                           /**< Mass 2 in solar masses */
                                          REAL8 r,                            /**< distance in metres*/
                                          LALDict *LALparams                     /**< LAL dictionary containing accessory parameters */
                                          )
{
    // lookup eco params
    REAL8 f0_rmodes1 = XLALSimInspiralWaveformParamsLookupECOFrequency1(extraParams);
    REAL8 delta_phi_rmodes1 = XLALSimInspiralWaveformParamsLookupECOPhaseDiff1(extraParams);
    REAL8 f0_rmodes2 = XLALSimInspiralWaveformParamsLookupECOFrequency2(extraParams);
    REAL8 delta_phi_rmodes2 = XLALSimInspiralWaveformParamsLookupECOPhaseDiff2(extraParams);
    
    REAL8 f0, f, df;
    UINT4 len, i;
    len = (*hptilde)->data->length;
    f0 = (*hptilde)->f0;
    df = (*hptilde)->deltaF;
    
    for (i=0; i<len; i++) { // loop over frequency points in sequence
        f = f0 + i*df;
        REAL8 rmode_phase_corr = 0.;
        COMPLEX16 rmode_phase_prefactor = 1. + 0.*I;

        if(f>=f0_rmodes1 && f<f0_rmodes2){
             rmode_phase_corr = (f/f0_rmodes1 - 1.)*delta_phi_rmodes1;
        }
        else if(f>=f0_rmodes2 && f<f0_rmodes1){
             rmode_phase_corr = (f/f0_rmodes2 - 1.)*delta_phi_rmodes2;
        }
        else if(f>=f0_rmodes1 && f>=f0_rmodes2){
             rmode_phase_corr = (f/f0_rmodes1 - 1.)*delta_phi_rmodes1 + (f/f0_rmodes2 - 1.)*delta_phi_rmodes2;
        }
        else {
             rmode_phase_corr = 0.;
        }
        rmode_phase_prefactor = cos(rmode_phase_corr) + I*sin(rmode_phase_corr);

        ((*hptilde)->data->data)[i] *= rmode_phase_prefactor;
        ((*hctilde)->data->data)[i] *= rmode_phase_prefactor;
    }

  return XLAL_SUCCESS;
}

```

1 place in LALSimInspiral.h:
```c
// around line ~920
/* ECO resonance */
int XLALSimECOResonanceTerm(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, REAL8 m1, REAL8 m2, REAL8 r, LALDict *LALparams);
```

3 places in LALSimInspiralWaveformParams.c: 
```c
DEFINE_INSERT_FUNC(EnableECO, INT4, "ECO", 0)
DEFINE_INSERT_FUNC(ECOFrequency1, REAL8, "ECOFreq1", 0)
DEFINE_INSERT_FUNC(ECOPhaseDiff1, REAL8, "ECOPhaseDiff1", 0)
DEFINE_INSERT_FUNC(ECOFrequency2, REAL8, "ECOFreq2", 0)
DEFINE_INSERT_FUNC(ECOPhaseDiff2, REAL8, "ECOPhaseDiff2", 0)

// the same for DEFINE_LOOKUP_FUNC and DEFINE_ISDEFAULT_FUNC
// lines 138, 280, 413

```

3 places in LALSimInspiralWaveformParams.h: 
```c
// around line ~125 (old v line ~27)
int XLALSimInspiralWaveformParamsInsertEnableECO(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertECOFrequency1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertECOPhaseDiff1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertECOFrequency2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertECOPhaseDiff2(LALDict *params, REAL8 value);
// around line ~253 (old v line ~110)
INT4 XLALSimInspiralWaveformParamsLookupEnableECO(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupECOFrequency1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupECOPhaseDiff1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupECOFrequency2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupECOPhaseDiff2(LALDict *params);
// around line ~380 (old v line ~200)
int XLALSimInspiralWaveformParamsEnableECOIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsECOFrequency1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsECOPhaseDiff1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsECOFrequency2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsECOPhaseDiff2IsDefault(LALDict *params);
```

## Checking ECO Changes in LALSimulation

```python
import lalsimulation as lalsim
import lal

extraParam = lal.CreateDict()
lalsim.SimInspiralWaveformParamsInsertEnableECO(extraParam, 1)
lalsim.SimInspiralWaveformParamsInsertECOFrequency1(extraParam, 128) #setting resonance frequency
lalsim.SimInspiralWaveformParamsInsertECOPhaseDiff1(extraParam, 100) #setting phase shift

#check the values are saved with Lookup functions
#useful in ipython
lalsim.SimInspiralWaveformParamsLookupECOFrequency1(extraParam)
lalsim.SimInspiralWaveformParamsLookupECOPhaseDiff1(extraParam)

#set variables
m1 = 28*lal.MSUN_SI
m2 = 23*lal.MSUN_SI
f_ref = 100.
phiRef = 0.
incl = 0.
s1x = 0.
s1y = 0.
s1z = 0.
s2x = 0.
s2y = 0.
s2z = 0.

phic = 0.
deltaF = 1./4.
f_min = 20.
f_max = 1024.
distance = 100e6*lal.PC_SI

#PhenomP arguments

model_flag = lalsim.IMRPhenomPv2_V

chi1_l, chi2_l, chip, thetaJN, alpha0, phi_aligned, zeta_polariz = lalsim.SimIMRPhenomPCalculateModelParametersFromSourceFrame(
                                                                   m1, m2, f_ref, phiRef, incl, s1x, s1y, s1z, s2x, s2y, s2z, lalsim.IMRPhenomPv2_V)

output_diff = lalsim.SimIMRPhenomP(chi1_l, chi2_l, chip, thetaJN, m1, m2, distance, 
                                  alpha0, phic, deltaF, f_min, f_max, f_ref, model_flag, extraParam_diff)
```