# Generating CompOSE HS(DD2) EOS Tables

This guide gives a complete, practical workflow for generating CompOSE EOS tables for the **HS(DD2) with electrons** EOS.

It covers two outputs:

1. a full 3D HDF5 table in `(T, n_b, Y_q)`, and
2. a cold beta-equilibrium slice.

The same workflow can be adapted to other CompOSE EOS tables by changing the EOS download link and using the ranges from that EOS’s `eos.t`, `eos.nb`, and `eos.yq` files.

---

## 1. What you will produce

After the full 3D workflow:

```text
eos.table.HSDD2.full3d.h5
```

After the beta-equilibrium workflow:

```text
eos.table.HSDD2.beta_T0p1MeV.dat
```

Optional one-column beta-equilibrium charge/electron fraction file:

```text
Ye_beq_HSDD2.txt
```

---

## 2. Download the CompOSE code

Create a working directory:

```bash
mkdir -p ~/compose_work
cd ~/compose_work
```

Clone the CompOSE source code:

```bash
git clone https://gitlab.obspm.fr/data_and_software_compose/code-compose.git
cd code-compose
```

If cloning with `git` is not available, download the source archive instead:

```bash
mkdir -p ~/compose_work
cd ~/compose_work

curl -L -o code-compose-master.tar.gz \
  https://gitlab.obspm.fr/data_and_software_compose/code-compose/-/archive/master/code-compose-master.tar.gz

tar -xzf code-compose-master.tar.gz
cd code-compose-master
```

In the commands below, replace `~/compose_work/code-compose` with your actual CompOSE directory if it is different.

---

## 3. Download the HS(DD2) EOS files

The HS(DD2) CompOSE page is:

```text
https://compose.obspm.fr/eos/18
```

Download the EOS archive into the CompOSE source directory:

```bash
cd ~/compose_work/code-compose

wget -O hsdd2_eos.zip \
  'https://compose.obspm.fr/download//3D/Hempel_SchaffnerBielich/hs_dd2_compose/with_electrons/eos.zip'
```

Unpack it:

```bash
unzip hsdd2_eos.zip
```

You should now have files like:

```text
eos.thermo
eos.compo
eos.micro
eos.t
eos.nb
eos.yq
eos.pdf
```

Check:

```bash
ls -lh eos.*
```

---

## 4. Check the EOS grid ranges

The files `eos.t`, `eos.nb`, and `eos.yq` contain the table grid. Use this command to print the number of points and the min/max values:

```bash
for f in eos.t eos.nb eos.yq; do
  awk 'NR==2{n=$1} NR==3{min=$1} NR>2{max=$1}
       END{print FILENAME, "npoints =", n, "min =", min, "max =", max}' "$f"
done
```

For HS(DD2), you should get:

```text
eos.t   npoints = 81   min = 1.0000000E-01   max = 1.5848932E+02
eos.nb  npoints = 326  min = 1.0000000E-12   max = 1.0000000E+01
eos.yq  npoints = 60   min = 1.0000000E-02   max = 6.0000000E-01
```

The first two lines of these files are metadata. The second line is the number of grid points, not a physical value.

---

## 5. Build CompOSE

### 5.1 Build with HDF5 support

The full 3D table should be generated as HDF5. On a cluster, load an HDF5 module first. The module name depends on the machine. For example:

```bash
module load cray-hdf5
```

Then check that the HDF5 Fortran compiler wrapper is available:

```bash
which h5fc
```

Build CompOSE with HDF5 enabled:

```bash
cd ~/compose_work/code-compose

unset USE_HDF5
make clean
rm -f *.o *.mod compose

make USE_HDF5=1 FC=h5fc compose 2>&1 | tee build_hdf5.log
```

Check that HDF5 was actually used:

```bash
grep -E "h5fc|have_hdf5|modhdf5|hdf5compose" build_hdf5.log
```

You want to see lines containing:

```text
h5fc
-Dhave_hdf5
modhdf5.f90
hdf5compose.f90
```

If the build uses only `gfortran` and does not mention `modhdf5.f90` or `hdf5compose.f90`, then HDF5 was not enabled.

### 5.2 Build without HDF5

For ASCII-only beta-equilibrium tables, a simple build can also work:

```bash
make clean
make compose
```

But for the full 3D table used by an HDF5 reader, use the HDF5 build above.

---

## 6. Generate the full 3D HDF5 table

The full 3D table depends on:

```text
T, n_b, Y_q
```

This is the table used by a 3D tabulated EOS reader.

---

### 6.1 Generate `eos.quantities`

Run:

```bash
./compose
```

Choose task:

```text
1
```

If it asks whether to generate a new `eos.quantities`, enter:

```text
1
```

CompOSE expects numbers at prompts. Do not type `y` or `yes`.

#### Regular thermodynamic quantities

When asked:

```text
How many regular thermodynamic quantities do you want to select for the file eos.table?
```

enter:

```text
7
```

Then enter these indices, one at a time:

```text
1
2
3
4
5
7
12
```

They mean:

| Index | Quantity |
|---:|---|
| 1 | pressure `p` |
| 2 | entropy per baryon `S` |
| 3 | shifted baryon chemical potential `mu_b - m_n` |
| 4 | charge chemical potential `mu_q` |
| 5 | lepton chemical potential `mu_l` |
| 7 | scaled internal energy per baryon `E/m_n - 1` |
| 12 | sound speed squared `c_s^2` |

#### Free-energy derivatives

Enter:

```text
0
```

#### Composition particles

CompOSE lists particle numbers and particle indices. For HS(DD2), the relevant entries are:

| Particle number | Particle | Particle index |
|---:|---|---:|
| 2 | neutron | 10 |
| 3 | proton | 11 |
| 4 | alpha | 4002 |

When asked:

```text
How many particles do you want to select for the file eos.table?
```

enter:

```text
3
```

Then CompOSE asks for particle numbers. Enter:

```text
2
3
4
```

Do not enter `10`, `11`, and `4002` here. Those are particle indices, not particle numbers.

The electron particle can be selected if you need it for diagnostics, but it is not needed for the usual reader setup because `Y_q` is already a table coordinate.

#### Average nucleus set

When asked:

```text
How many sets do you want to select for the file eos.table?
```

enter:

```text
1
```

If CompOSE asks for the set number, enter:

```text
1
```

This selects the average nucleus set with index `999`.

#### Microscopic quantities

Enter:

```text
0
```

#### Error estimates

Enter:

```text
0
```

#### Output format

If CompOSE asks:

```text
1: ASCII, else: HDF5
```

enter:

```text
2
```

This selects HDF5 output.

---

### 6.2 Generate `eos.parameters`

Run:

```bash
./compose
```

Choose task:

```text
2
```

#### Interpolation order

When asked for the interpolation order in `T`, `n_b`, and `Y_q`, enter:

```text
1
1
1
```

#### Beta equilibrium

When asked:

```text
Please select if you want to calculate the EoS of matter in beta-equilibrium.
1: yes, else: no
```

enter:

```text
0
```

#### Fixed entropy

When asked:

```text
Please select if you want to calculate the EoS for given entropy per baryon.
1: yes, else: no
```

enter:

```text
0
```

#### Tabulation scheme

When asked:

```text
Please select the tabulation scheme for the parameters from
0: explicit listing of parameter values
1: loop form of parameter values
```

enter:

```text
1
```

#### Temperature grid

Enter:

```text
1.0000000E-01 1.5848932E+02
81
1
```

This means:

```text
T_min T_max = 1.0000000E-01 1.5848932E+02
number of T points = 81
logarithmic scaling = 1
```

#### Baryon density grid

Enter:

```text
1.0000000E-12 1.0000000E+01
326
1
```

This means:

```text
n_b,min n_b,max = 1.0000000E-12 1.0000000E+01
number of n_b points = 326
logarithmic scaling = 1
```

#### Charge fraction grid

Enter:

```text
1.0000000E-02 6.0000000E-01
60
0
```

This means:

```text
Y_q,min Y_q,max = 1.0000000E-02 6.0000000E-01
number of Y_q points = 60
linear scaling = 0
```

CompOSE should now write:

```text
eos.parameters
```

---

### 6.3 Generate the table

Run:

```bash
./compose
```

Choose task:

```text
3
```

A successful HDF5 run should include messages like:

```text
format of output table = 2 (HDF5)
call writing HDF5 table
writing 7 thermodynamic quantities into file
writing 3 pairs into file
writing 1 quadruples into file
```

For HDF5 output, CompOSE usually writes the actual HDF5 file as:

```text
eoscompose.h5
```

Copy it to a clearer name:

```bash
cp eoscompose.h5 eos.table.HSDD2.full3d.h5
```

---

### 6.4 Check the HDF5 file

Check that it is really HDF5:

```bash
file eos.table.HSDD2.full3d.h5
```

Expected:

```text
Hierarchical Data Format (version 5) data
```

List the groups:

```bash
h5ls eos.table.HSDD2.full3d.h5
```

Expected groups include:

```text
Composition_pairs
Composition_quadrupels
Parameters
Thermo_qty
metadata
```

Check the selected quantities:

```bash
h5dump -d /Thermo_qty/index_thermo eos.table.HSDD2.full3d.h5
h5dump -d /Composition_pairs/index_yi eos.table.HSDD2.full3d.h5
h5dump -d /Composition_quadrupels/index_av eos.table.HSDD2.full3d.h5
```

Expected values:

```text
/Thermo_qty/index_thermo         = 1, 2, 3, 4, 5, 7, 12
/Composition_pairs/index_yi      = 10, 11, 4002
/Composition_quadrupels/index_av = 999
```

---

## 7. Generate a cold beta-equilibrium slice

The beta-equilibrium slice is a 1D table along `n_b`. It is useful for cold EOS checks and TOV-style preprocessing.

For HS(DD2), use the lowest table temperature:

```text
T = 1.0000000E-01 MeV
```

Before changing the CompOSE input files, save the full 3D setup:

```bash
cp eos.quantities eos.quantities.full3d
cp eos.parameters eos.parameters.full3d
```

---

### 7.1 Generate beta-slice `eos.quantities`

Run:

```bash
./compose
```

Choose task:

```text
1
```

Regenerate `eos.quantities`:

```text
1
```

For a minimal beta-equilibrium table, select three regular thermodynamic quantities:

```text
3
21
1
3
```

These are:

| Index | Quantity |
|---:|---|
| 21 | internal energy per baryon `E` [MeV] |
| 1 | pressure `p` |
| 3 | shifted baryon chemical potential `mu_b - m_n` |

For the remaining groups, enter:

```text
0
0
0
0
0
```

That means:

```text
free-energy derivatives: 0
composition particles:   0
average sets:            0
microscopic quantities:  0
error estimates:         0
```

For a simple beta-equilibrium text table, choose ASCII output:

```text
1
```

If you want extra diagnostics, choose six regular quantities instead:

```text
6
21
1
3
12
15
2
```

This gives:

```text
E, p, mu_b - m_n, c_s^2, Gamma, entropy
```

---

### 7.2 Generate beta-slice `eos.parameters`

Run:

```bash
./compose
```

Choose task:

```text
2
```

Use interpolation order 1 in all variables:

```text
1
1
1
```

Turn on beta equilibrium:

```text
1
```

Turn off fixed entropy:

```text
0
```

Use loop form:

```text
1
```

Use one temperature point at the lowest table temperature:

```text
1.0000000E-01 1.0000000E-01
1
0
```

Use the full baryon density range:

```text
1.0000000E-12 1.0000000E+01
326
1
```

Use the full charge-fraction range as the search range for the beta-equilibrium solve:

```text
1.0000000E-02 6.0000000E-01
```

If CompOSE asks for a number of `Y_q` points and scaling, enter:

```text
60
0
```

If it does not ask for those two values in beta-equilibrium mode, that is fine.

---

### 7.3 Generate the beta-equilibrium table

Run:

```bash
./compose
```

Choose task:

```text
3
```

Rename the output immediately:

```bash
mv eos.table eos.table.HSDD2.beta_T0p1MeV.dat
cp eos.beta eos.beta.HSDD2_T0p1MeV.dat 2>/dev/null || true
cp eos.report eos.report.HSDD2.beta_T0p1MeV.txt 2>/dev/null || true
```

Check:

```bash
head eos.table.HSDD2.beta_T0p1MeV.dat
wc -l eos.table.HSDD2.beta_T0p1MeV.dat
```

For the minimal quantity selection, the columns are:

```text
1  T [MeV]
2  n_b [fm^-3]
3  Y_q
4  E [MeV]
5  p [MeV fm^-3]
6  mu_b - m_n [MeV]
```

---

## 8. Create a one-column beta-equilibrium `Y_e` file

For charge-neutral matter in this table, the beta-equilibrium `Y_e`/charge fraction is the third column of the beta table.

Create a one-column file:

```bash
awk '!/^#/ {printf "%.18e\n", $3}' \
  eos.table.HSDD2.beta_T0p1MeV.dat > Ye_beq_HSDD2.txt
```

Check it:

```bash
head Ye_beq_HSDD2.txt
wc -l Ye_beq_HSDD2.txt
```

For the full HS(DD2) `n_b` grid, the file should have:

```text
326
```

lines.

---

## 9. Optional: extract a beta slice from the full 3D HDF5 table

If you already have the full 3D HDF5 table, you can also post-process it by solving:

```text
mu_l(T, n_b, Y_q) = 0
```

along the `Y_q` direction.

If you have the helper script `generate_beta_eq_slice.py`, run:

```bash
python3 generate_beta_eq_slice.py eos.table.HSDD2.full3d.h5 \
  --cold \
  --out-prefix HSDD2_beta_cold
```

Then make the one-column file:

```bash
awk '!/^#/ {printf "%.18e\n", $3}' HSDD2_beta_cold.dat > Ye_beq_HSDD2.txt
```

This method is useful for quick post-processing. The built-in CompOSE beta-equilibrium mode is the preferred method for production beta slices.

---

## 10. Restore the full 3D setup

If you saved the full 3D input files and want to regenerate the full 3D table later:

```bash
cp eos.quantities.full3d eos.quantities
cp eos.parameters.full3d eos.parameters

./compose
# choose task 3
```

---

## 11. Notes for EOS readers

A typical HDF5 reader should use:

```text
/Parameters
/Thermo_qty
/Composition_pairs
/Composition_quadrupels
```

The full 3D table generated above contains:

```text
index_thermo = 1, 2, 3, 4, 5, 7, 12
index_yi     = 10, 11, 4002
index_av     = 999
```

For chemical potentials:

```text
3 = mu_b - m_n
4 = mu_q
5 = mu_l
```

A common mapping is:

```text
mu_n = mu_b - m_n
mu_p = mu_n + mu_q
mu_e = mu_l - mu_q
```

If your reader has a variable for `mu_l` or `MUNU`, fill it from CompOSE index `5`.

---

## 12. Quick troubleshooting

### The generated `.h5` file is not HDF5

Check:

```bash
file eos.table.HSDD2.full3d.h5
```

If it says `ASCII text`, the file is ASCII even if it has a `.h5` name. Regenerate with HDF5 output selected.

A real HDF5 file should say:

```text
Hierarchical Data Format (version 5) data
```

### HDF5 output says it wrote `eos.table`, but there is no `eos.table`

For HDF5 output, the actual file is usually:

```text
eoscompose.h5
```

Use:

```bash
ls -lh eoscompose.h5
cp eoscompose.h5 eos.table.HSDD2.full3d.h5
```

### The build did not enable HDF5

Check the build log:

```bash
grep -E "h5fc|have_hdf5|modhdf5|hdf5compose" build_hdf5.log
```

If those strings are missing, rebuild with:

```bash
make clean
rm -f *.o *.mod compose
make USE_HDF5=1 FC=h5fc compose 2>&1 | tee build_hdf5.log
```

### CompOSE crashes after typing `yes`

CompOSE prompts expect numbers, not words. Use:

```text
1
```

instead of:

```text
yes
```

### Files get overwritten

CompOSE reuses filenames such as:

```text
eos.quantities
eos.parameters
eos.table
eoscompose.h5
eos.report
```

Save important files before switching workflows:

```bash
cp eos.quantities eos.quantities.full3d
cp eos.parameters eos.parameters.full3d
cp eoscompose.h5 eos.table.HSDD2.full3d.h5
```

---

## 13. Final checklist

For the full 3D table:

```bash
file eos.table.HSDD2.full3d.h5
h5ls eos.table.HSDD2.full3d.h5
h5dump -d /Thermo_qty/index_thermo eos.table.HSDD2.full3d.h5
h5dump -d /Composition_pairs/index_yi eos.table.HSDD2.full3d.h5
h5dump -d /Composition_quadrupels/index_av eos.table.HSDD2.full3d.h5
```

Expected:

```text
index_thermo = 1, 2, 3, 4, 5, 7, 12
index_yi     = 10, 11, 4002
index_av     = 999
```

For the beta-equilibrium slice:

```bash
head eos.table.HSDD2.beta_T0p1MeV.dat
awk '!/^#/ {printf "%.18e\n", $3}' \
  eos.table.HSDD2.beta_T0p1MeV.dat > Ye_beq_HSDD2.txt
wc -l Ye_beq_HSDD2.txt
```
