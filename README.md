# Accurate pseudopotentials for metagga calculations
**D. R. Hamann**

_Department of Physics and Astronomy, Rutgers University, Piscataway, New Jersey and
Mat-Sim Research LLC, Lafayette, Colorado_

May 16, 2024

## Introduction
METAPSP-1.0.1 is a specialized adaptation of ONCVPSP-4.0.1[<sup>1</sup>](#ref-1) developed to deal with metagga functionals. It has proven possible to generate pseudopotentials satisfying the conditions of generlized norm conservation[<sup>2</sup>](#ref-2) despite the more complicated Kohn-Sham equations required to minimize these functionals. In particular, M101 atomic pseudo wavefunctions exactly reproduce the corresponding all-electron wavefunctions outside their "core radii" $r_c$. Their norms (integrated charge densities inside $r_c$), their radial log derivatives at $r_c$, as well as the first energy derivatives of these log derivatives match the all-electron values at atomic bound states energies and those of selected scattering-states. The convergence of the Fourier transforms of these pseudo wave functions is optimized in the same fashion as those of their lda-gga predecessors in ON401.[<sup>3</sup>](#ref-3) These meta psps will, however, only be of value with application codes that fully treat the Kohn-Sham complications.

The code for performing the mgga calculations is sufficiently different from that implemented for lda/gga psps in ON401 that new routines have been introduced to replace many of those. Most have predecessors *.f90 from which they are distinguished as *_m.f90. The key differences are discussed below. For a preview of the nice pseudopotentials M101 can produce, see the graphics showing a W metapsp treating 8 states as valence in [W_metapsp_plots.pdf](W_metapsp_plots.pdf).

The new code can be built in the manner familiar from ON401. The active library links in [make.inc](make.inc) must be edited to suit local libraries, etc.,and the command make onvpspm issued in src will build the executable oncvpspm.x. The simple shell scripts in the directory test will execute the calculations for a specified file of input data, and automatically generate diagnostic plots using gnuplot.

**Radial Equations**

Because of the dependence of the mgga functionals on the kinetic energy density $\tau(\mathbf{r})$ the Kohn-Sham equations contain a gradient term multiplied by a vector potential $\mathbf{v}_\tau(\mathbf{r})$. While independently derived, our radial Schroedinger equation for a spherical atom is identical to Eqs. 21 & 22 of Holzwarth et al.[<sup>4</sup>](#ref-4) and contains a scalar function denoted vtau multiplying a radial first derivative. These terms are used for the all-electron atom and for the pseudo atom.

**MGGA functional**

Despite my initial intention of linking to libxc for a choice of funtionals, extensions of the ON401 interface code retained unresolved errors. I chose to proceed with Pendry's R2SCAN functional,[<sup>5</sup>](#ref-5) for which Natalie Holtzwarth generously provided me local routines that could be directly incorporated into the M101 source directory. With this in hand, Natalie, Susi Letola and I cross-checked our calculations of several atoms and found agreement in eigenvalues and total energies to within numerical accuracy. Three different formulations were used for solution of the radial equations, so I'm confident of the accuracy of the code presented here.

I will be very appreciative of users of the present code who could contribute a libxc interface which can reproduce the local R2SCAN01 results, and hence presumably those of other functionals. On the other hand, since even this well-behaved and internally smoothed functional presented various convergence stability problems, only a limited group of others might be usable with the present algorithms. I've been persuaded that R2SCAN01 is a very sophisticated functional, and testing to assess the value of these M101 psps in applications programs is a more important next step.

**All-electron atom**

The overall calculation is initialized by first performing a PBE functional gga calculation for the atom to obtain a good starting potential and trial set of eigenvalues for the mgga atom. Its vtau is initialized as zero. Despite the good start, the Anderson iteration scheme used for lda/gga proved unstable, and was replaced by simple mixing with a small parameter.

The algorithm used to solve the radial equation in the routine [lschfb_m](src/lschfb_m.f90) is similar to that used in oncvpsp. The trial eigenvalue is iterated until the log derivatives of outward and inward integrations of the wave function match at the classical turning point. A logarithmic radial grid is utilized. The second derivative of the wave function at each new mesh point is calculated from the extrapolated values of the first and second derivatives multiplied by appropriate coefficients, and the extrapolations are corrected by several iterations at that point. Including vtau times the gradient adds a term to each of the coefficients multiplying those first and second derivatives. This algorithm is stable, however the peculiar behavior of the r2scan01 exchange potential at large radii led to serious errors in identifying the classical turning point, and the search for this point had to be modified.

Scalar relativistic corrections to the radial equation^6 also make additions to those two coefficients. As in ON401, exact generalized norm conservation can no longer be preserved, but users have been happy to tolerate the resulting ~$10^{-4}$ errors. In the meta case, these additional terms must be smoothly cut off approaching the smallest core radius, so there is already an error, usually of even smaller magnitude, in the all-electron reference results. As far as I know, relativistic forms of the mgga total energy do not exist, and there is no systematic derivation of such corrections. So, the same additive terms used in lda/gga calculation were simply added. For the time being, spin-orbit terms calculated with ON401 using PBE or some other gga functional should probably provide reasonable accuracy for fully-relativistic mgga psp calculations.

At this point the real trouble started. The all-electron atoms no longer converged, and noise on the log-mesh-spacing scale was found in the first few decades of the mesh. Its amplitude grew exponentilly with iterations. I speculate that this trouble was caused by the fact that, effectively, a third derivative of a function of the density was making its was into the above-mentioned coefficients through the r2scan01 exchange potential. My ad-hoc solution was to smooth that exchange potential (and vtau for good measure) by Gaussian convolution over the first few decades, smoothly crossing over the the originals at larger radii. The Gaussian width was a specified number of mesh points, typically 25, with the Gaussian equal to $10^{-2}$ at the first and last points. The effects of this smoothing were checked against unsmoothed results for non-relativistic calculations, and were identical within numerical accuracy for eigenvalues and total energies.

**All-electron "scattering" states**

Pseudopotentials satisfying GNC at positive energies were generated from bound-state all-
electron wave functons in specially-constructed "confining well" potentials which were exactly equal to the full potential and vtau inside $r_c$ for each $\ell$. While the corresponding wellstate.f90 routine in the ON401 code simply added a potential that went to a fixed constant value above the desired energy, the often odd behavior of the r2scan01 vxc and vtau in the atomic tails ultimately led me to smoothly kill the AE vxc and vtau outside the largest $r_c$, while adding the smoothly-rising barrier. The algorithms for adjusting the radial scale and amplitude of that barrier in wellstate_m.f90 ended up being considerably more complex than those in the ON401 routine.

**Pseudo-wave-functions**

The pswfs are constructed essentially exactly as in ON401. They are expressed as linear
combinations of spherical Bessel functions inside $r_c$. Their number is adequate to match the value, the specified number of derivatives, and the integrated charge density of the corresponding AE funcion, to provide orthogonality to the lower-energy pswfs and to still have a few left over to optimize the pswf Fourier convergence. The largest wavevector of the sbfs is qcut specified in the data. Outside $r_c$ they are identical to the AE functions.

One feature has been added. Optimizing the Fourier convergence essentially means moving as much of the pseudo charge density as the constraints allow outwards towards $r_c$. There are no constraints in this algorithm to keep the first pswf above zero near the origin, and and as users of ON401 know, the message "ERROR first l=? pseudo wave function has node" required input data changes. This condition must be satisfied to construct the pseudopotential. The most common input change is to reduce qcut. So many such errors occurred using my large test-set of existing input data files that the meta code now automatically reduces qcut by 5% increments six times. If the node still
occurs, then an ERROR exit is generated, and manual reduction of the input $r_c$ is usually the next step. M101 users will appreciate this automation, which eliminated almost all the ERROR exits.

**Local Potential and "Core" models**

In ON401, the local potential accompanying the non-local projectors in the complete
pseudopotential could be chosen as the smooth analytic continuation of the AE potential to the origin, or the so-called "semi-local" potential for a chosen $\ell$. The latter choice essentially meant using the original norm-conserving psp[<sup>7</sup>](#ref-7) and forgoing accuracy at any higher bound- or scattering-state energies for this $\ell$. This choice is no longer available for M101 because similarly-constructed functions do not account for meta issues and are discontinuous at $r_c$. So the input data choice lloc=4, which indicates the analytic continuation approach, is the only M101 option. Despite the discontinuities, these functions are plotted as "Unscreened Semi-Local Guides" in the automated graphics output, because a good choice of input data parameters lpopt, rc(5), and dvloc0 should put the local potential close them.

Mgga also requires a $\mathbf{v}_\tau$ to be used in the pseudopotential which is identical to the AE $\mathbf{v}_\tau$ from the smallest $r_c$ out. The output local potential is the unscreened local potential, which is "re-screened" when used in applications by the Hartree and exchange-correlation potentials calculated there. There is no such analog for an "unscreened local $\mathbf{v}_\tau$," since $\mathbf{v}_\tau$ comes strictly from the mgga functional. Instead, the idea of the nonlinear core correction[<sup>8</sup>](#ref-8) is extended in M101. Whereas the ON401 code jumps directly from pseudo wave function optimization to non-local projector construction (run_vkb.f90), M101 first models $\rho$s and $\tau$s which, when added to the pseudo $\rho$ and $\tau$, will reproduce the AE rho and tau from the minimum rc out. Outside the maximum rc these will match the AE core values, so this construction has some relation to the older approach.

As discussed previously, convergence optimization shifts the weight of pseuo wave functions towards $r_c$. The AE – PS rho and tau differences usually become negative for some ways inside the the maximum $r_c$ until the core rho and tau take over. Continuity requires models to follow this negative behavior a limited distance inside the minimum $r_c$. The most effective strategy I've developed proceeds in two steps, the same for both rho and tau, and is specified by icmod=5 in the input data. The maximun values of the PS rho and tau are determined. A polynomial with the values of the datum rcfact times those maxima at $r=0$ and going to zero with 3 zero derivatives at rcmin is added to the differences. Next, a search is made for that radius inside rcmin from which an analytic continuation of the above functions matching their values and 3 derivatives can be made to achieve values of fcfact times the respective maxima at the origin. This sounds quite complicated, but the choices rcfact=0.5 and fcfact=3, used by default unless icmod already equals 5 in the input data, typically produce "core" models very similar to model rhos adjusted by hand for ON401. It is highly improbable that the small negative regions that often occur in the new models near rcmin will result in net negative values once added to the pseudo $\rho$ and $\tau$ in any application.

**Non-local Projectors**

It is the next step in M101 that contains the significant advance making these calculations of value: a formulation for the non-local projectors that incorporates the effects of the vtau term in the radial equation in such a way that generalized norm conservation can be achieved.[<sup>9</sup>](#ref-1) That this is the correct formulation has only been established by what I'll call numerical experiments, not by any analytic derivation. The "Diagnostics" output and log-derivative results achieve a level of accuracy limited only by numerical convergence. The spacing of the logarithmic radial mesh was reduced by a factor of 4 to obtain convergence comparable to the lda/gga calculations. As in ON401, the inclusion of scalar-relativistic calculations continues to introduce errors of order $10^{-4}$. Unlike the meta changes to the radial equations, these are not carried over to pseudopotential calculations, and thus cannot be corrected by further changes to the projectors.

**Diagnostics, Ghost detection, etc.**

The [diagnostics_m.f90](src/run_diag_m.f90) routine is a straightforward generalization of that in ON401, simply substituting the meta-adapted all-electron and pseudo radial equations.

"Ghosts" are unintended solutions of the pseudopotential radial equations. The non-local
projectors permit solutions violating the usual ordering of energies with node count. They are typically highly localized around those projectors, and are in no way unique to the metagga complications. The routine [run_ghosts_m.f90](src/run_ghosts_m.f90) very much follows the ON401 approach. The screened local pseudopotential is terminated by a hard wall at $3r_c$. The meta radial equation with the pseudo vtau but no non-local
projectors is solved for a set of basis functions ranging up to an energy of approximately 1.25 times the the kinetic energy of the highest-wavevector spherical Bessel function. The non-local potential is then diagonalized in that basis. Any solutions at energies lower than those corresponding to the atomic bound states are identified as negative-energy ghosts, and a WARNING is printed in the output.

Any positive-energy solutions whose mean radii are less than rc are likely to give rise to spurious very narrow conduction bands, and a GHOST warning is issued. This identification is less certain, and the automatically-generated log-derivative plots comparing all-electron and pseudopotential results should be examined up to the highest energy of interest. The confined AE potential might also have such localized states within the hard wall, but not create spurious bands in an application. Such a report of ghosts which happen to be within the range of the specified log-derivative plots and those plots are shown in the [Ghosts.pdf](Ghosts.pdf) file.

The routines [run_plot_m.f90](src/run_plot_m.f90), [run_pshft_m.f90](src/run_pshft_m.f90), and [gnu_script_m.f90](src/gnu_script_m.f90) are straightforward generalizations of the ON401 versions. Two new plots have been added. The all-electron and model plus pseudo taus are shown with a log vertical scale following the rho plot. Following that, the AE and M+PS vtaus are compared. A log horizontal scale is used for these to illustrate the potentially large effects of vtau in low-density regions of materials. It will be interesting to see if this behavior is particular to R2SCAN01 or arises in other meta psps.

**Application output**

As all the above discussion indicates, the application needs to accept the model $\tau$ input as well as the model $\rho$. The output $\rho_{ps}$, which can be used for initialization is in the ON401 output, is supplemented by the $\tau_{ps}$ in M101. The local potential as well as the non-local projectors remain as in ON401, as do the radial grid output and various "bookkeeping" data.

The last entry on the first active line of the input data is psfile, indicating the application output format. The accepted values for this text variable are 'none', 'psp8', 'upf', and 'psml'. The 'none' entry is useful for developing and tuning data when all you really want to see is the diagnostics and plots.

The psp8 data format, originally developed for Abinit, is unchanged in structure. However, in M101 the last two columns of text previously used for the 3rd ahd 4th derivatives of the model core rho contain the pseudo and model taus, respectively. The higher derivatives of the model core rho can be generated internally by Abinit when higher-order density-functional perturbation theory is implemented for meta psps (don't hold your breath!). The upf format initially developed for Quantum Espresso has been extended to include the ps and model taus with the gracious help of Paolo Giannozzi, Matthieu Verstraete, and Pietro Delugas. These are in separate, appropriately labeled data blocks, and the fact that this is data for use with r2scan01 metagga calculations is indicated in header listings.

The psml data format was not included in my ON401 release, but was developed for the Siesta electronic structure code by Alberto Garcia and coworkers. He has helped me adapt their code for ON401 to include the ps and model taus plus appropriate header labeling. Building M101 requires first building and installing the xmlf90 library. I quote his instructions in the [xmlf90_build.txt](xmlf90_build.txt) file. Note that the radial grid for the psml output is not the uniformly-spaced grid of the psp8 and upf formats, but a quasi-uniform sampling of a fraction of the points of the original log grid.

If psml output is not needed, the command make oncvpspmx in the src directory will build the alternate executable oncvpspmx.x, and the libxml references in [make.inc](make.inc) can be omitted.


**Testing**

Nearing the final stages of development of M101, I assembled a collection of 380 input data
files developed for ON401 and earlier lda/gga releases. These comprised 91 from the SG
collection,[<sup>10</sup>](#ref-10) 178 from the Pseudo-Dojo collection,[<sup>11</sup>](#ref-11) and 111 from my own unpublished work (un-curated). Some of the latter were designed for high-temperature simulations, treating several core levels as valence, and ensuring accurate scattering up to 10s of Ha by employing 3 or 4 projectors per $\ell$. I had initiated the ONCVPSP project with the belief that those using psps should get some deeper understanding of them by making their own. However, the well-tested published psps far surpassed home-made ones in usage. The M101 code substituted defaults for input variables not consistent with its requirements. This entire set could be run with a high degree of parallelization in a few minutes. Scanning the outputs, truncated at various stages of development, provided feedback which led to error correction and various refinements of the code. Most of this development was done with non-relativistic calculations where exact agreement with AE was expected for a number of quantities.

For the final version released here, errors stopped the code in 2 cases because of first-pseudo-wavefunction nodes that were not removed by the automated reduction of the qcut input variable. One stopped with the unusual message "ERROR ps0norm > uunorm, code will stop." All these could be fixed by further changes in qcut, rc, or both. The scalar-relativistic runs using exactly the same data and defaults yielded only 1 "node" error. Lots of ghosts were found, and I'm not sure how those statistics would compare with lda/gga calculations or whether the default substitutions were at least in
part responsible. Many positive-energy ghosts can be ignored, and negative-energy ones eliminated, usually by tuning input parameters to somewhat relax localization and cutoff-convergence.

That so much existing data ran without errors did not necessarily mean that the resulting
r2scan01 psps were "good." Examination of a random selection of these outputs by re-plotting the automatic graphics revealed a number of cases where tuning the data could yield smoother non-local projectors, smaller negative regions in the model rhos and taus, etc. First-row atoms were particularly challenging in this regard.

As in ON401, a "tests" directory in this distribution contains sample input files in a "data" sub- directory, corresponding output files in "refs" and "refs_sr" directories. It includes simple shell scripts to test your build of the code by running the data. Input files from my own previous work were used for these, but are shown almost exclusively with the default substitutions copied from the output. Unlike ON401, M101 runs the same object code, oncvpspm.x, for both non-relativistic and scalar-
relativistic calculations, the latter option being chosen when the file SR (contents also "SR") is temporarily created in the running directory.

I'm releasing this code without ANY systematic testing of these psps' performance in
applications. Matthieu Verstraete has collaborated with me in developing an initial interface for Abinit and verifying that valance atomic levels could be reasonably reproduced for a few single atoms in large unit cells. Fully-converged testing on a number of cases lies in the future. More important, significantly improved agreement with experiment using M101 psps in codes running r2scan01 is not a-priori guaranteed. I've gathered that in general significant metagga improvements are found in AE solid calculations for weakly-bonded materials. Although I certainly hope not, it is plausible that gga psps could combine well enough with application metagga improvements in the treatment of low-density regions in solids to accomplish similar results.

<a name="ref-1"></a>
1. D. R. Hamann, Phys. Rev. B **88** , 085117 (2013).
<a name="ref-2"></a>
1. D. Vanderbilt, Phys. Rev. B **41** , 7892 (1990).
<a name="ref-3"></a>
1. A. M. Rappe _et al._ , Phys. Rev. B **41** , 1227 (1990).
<a name="ref-4"></a>
1. N. A. W. Holtzwarth _et al._ , Phys. Rev. B **105** , 125144 (2022).
<a name="ref-5"></a>
1. J. W. Furness _et al._ , Phys. Chem. Lett. **11** , 8208 (2020).
<a name="ref-6"></a>
1. D. D. Koeling and B. N. Harmon, J. Phys. C **10** , 3107 (1977).
<a name="ref-7"></a>
1. D. R. Hamann _et al._ , Phys. Rev. Lett. **43** , 1494 (1979).
<a name="ref-8"></a>
1. S. G. Louie _et al._ , Phys. Rev. B **43** , 1993 (1991).
<a name="ref-9"></a>
1. D. R. Hamann (unpublished).
<a name="ref-10"></a>
1.  M. Schlipf and F. Gygi, Computer Physics Communications **196** , 36 (2015).
<a name="ref-11"></a>
1.  M. J. van Setten _et al_ ., Computer Physics Communications **226** , 39 (2018)
