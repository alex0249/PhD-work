## Mechanics: stress and deformation

**Strain**
How much a material stretches, compresses, or changes shape relative to its original size.

For a change in length:
$$
\varepsilon = \frac{L-L_0}{L_0}
$$
$\varepsilon$ = engineering strain; $L_0$ = original length; $L$ = current length. Strain has no units.

**Stress**
Internal force per unit area within a material.

Average normal stress across an area:
$$
\sigma_{\mathrm{avg}} = \frac{F_n}{A}
$$
$F_n$ = force perpendicular to the area; $A$ = area carrying that force. Stress is measured in Pa. Local stress can vary across the area.

**Normal stress**
Stress acting perpendicular to a surface.

**Shear stress**
Stress acting along a surface.

**Hydrostatic stress**
The average of the three normal stresses; the part associated with uniform tension or compression.

$$
\sigma_h = \frac{\sigma_{xx}+\sigma_{yy}+\sigma_{zz}}{3}, \qquad p=-\sigma_h
$$
$\sigma_h$ = hydrostatic stress; $\sigma_{xx},\sigma_{yy},\sigma_{zz}$ = normal stresses in three perpendicular directions. Here tension is positive and pressure is positive in compression. Stack pressure is not generally equal to local $p$.

**Deviatoric stress**
The part of stress remaining after removing hydrostatic stress; associated with distortion.

$$
\mathbf{s}=\boldsymbol{\sigma}-\sigma_h\mathbf{I}
$$
$\mathbf{s}$ = deviatoric stress tensor; $\boldsymbol{\sigma}$ = total stress tensor; $\sigma_h$ = hydrostatic stress; $\mathbf{I}$ = identity tensor (ones on the diagonal, zeros elsewhere).

**Isochoric deformation**
Deformation that changes shape without changing volume.

**Dilatation**
Relative volume change: positive for expansion and negative for contraction. For small strains, it is the sum of the three normal strains.

**Densification**
An increase in density as the same amount of material occupies less volume. In glasses, the atomic network can pack more tightly without actual pores being present.

## Mechanics: material behaviour

![[Attachments/Concept diagrams/stress-strain-ductile.png|600]]

*Stress–strain curve:* the initial slope is Young's modulus; yielding marks significant plastic deformation. This is an illustrative engineering stress–strain curve for a ductile metal, not measured battery-material data. Brittle ceramics may fracture with little plastic deformation.

**Work hardening (strain hardening)** — Increased resistance to further plastic deformation as plastic strain increases, usually because dislocations multiply and obstruct one another. Annealing can limit this by rearranging or removing dislocations.
**Necking** — Local narrowing of a specimen during stretching.
**Ultimate tensile strength** — The maximum engineering tensile stress on the curve.

**Elasticity**
A material’s ability to recover its original shape and size when the force is removed.

**Elastic deformation**
The recoverable change in shape or size. For example, gently stretching a spring and releasing it.

**Plasticity**
A material’s ability to undergo permanent deformation without breaking.

**Plastic deformation**
The permanent change that remains after the force is removed. For example, bending a paperclip so it stays bent.

**Young's modulus**
A measure of stiffness: resistance to elastic stretching or compression.

For linear elastic loading in one direction, with the sides free:
$$
\sigma=E\varepsilon^e
$$
$E$ = Young's modulus (Pa); $\sigma$ = axial stress; $\varepsilon^e$ = elastic strain along the loading direction. Chemical or plastic strain is not included in $\varepsilon^e$.

**Poisson's ratio**
How much a material contracts sideways when stretched lengthwise.

For uniaxial elastic loading:
$$
\nu=-\frac{\varepsilon^e_{\mathrm{sideways}}}{\varepsilon^e_{\mathrm{lengthwise}}}
$$
$\nu$ = Poisson's ratio; the strains are the elastic changes perpendicular and parallel to the applied load. The minus sign accounts for sideways contraction during stretching.

**Isotropy**
Having the same material properties in every direction.

**Anisotropy**
Having material properties that depend on direction.

**Yielding**
The onset of significant plastic deformation.

**Viscoelastic strain**
Deformation combining elastic recovery with a time-dependent response. Unlike ideal elastic deformation, it depends on loading history and can recover gradually after unloading. Recovery may be complete, or the material/model may also include irreversible flow.

![[Attachments/Concept diagrams/elastic-viscous-viscoelastic-comparison.png|650]]

*For a load applied at $t_a$ and removed at $t_r$:*
- **Elastic:** strain appears immediately and disappears immediately after unloading.
- **Viscous:** strain accumulates while loaded and remains after unloading in this ideal model.
- **Viscoelastic:** combines immediate and gradual responses; recovery also takes time.

$t_a$ = time the load is applied; $t_r$ = time it is removed. The “+” is schematic: the exact curve depends on the material/model. Some viscoelastic materials recover completely; models including irreversible viscous flow retain strain.

Image supplied by you; original source not identified.

**Creep**
Deformation that grows over time under sustained stress. For example, lithium can slowly change shape under stack pressure, helping maintain contact with the electrolyte. Creep describes behaviour over time; it can arise from different mechanisms.

![[Attachments/Concept diagrams/creep-stages.jpg|450]]

*Top:* stress is applied and held constant. *Bottom:* strain jumps initially, then continues increasing through creep. The slope is the strain rate. $\sigma$ = applied stress; $\varepsilon$ = strain; $\varepsilon_0$ = initial strain just after loading; $t_0$ = loading time; $t$ = time.

**Primary creep** — Creep rate decreases with time.
**Secondary creep** — Creep rate is approximately constant.
**Tertiary creep** — Creep rate increases, often as damage develops towards failure.

This is a general schematic, not lithium data; not every material or battery condition shows all three stages.

A simple steady-creep model at a fixed temperature (Norton's law, using magnitudes in one-dimensional loading) is:
$$
\dot{\varepsilon}_{\mathrm{cr}}=\frac{d\varepsilon_{\mathrm{cr}}}{dt}=A|\sigma|^n
$$
$\varepsilon_{\mathrm{cr}}$ = creep strain; the dot = rate of change with time; $\dot{\varepsilon}_{\mathrm{cr}}$ = creep strain-rate magnitude (s$^{-1}$); $\sigma$ = applied axial stress (Pa); $A$ = fitted coefficient (Pa$^{-n}$ s$^{-1}$); $n$ = dimensionless stress exponent. $A$ and $n$ depend on the material and temperature. This model does not describe all creep stages; lithium parameters must match the relevant conditions.

**Stress relaxation**
Stress decreasing over time while deformation is held fixed.

**Maxwell relaxation time**
The characteristic time over which stress relaxes in a simple Maxwell viscoelastic model when total strain is held fixed. After one relaxation time, about 37% of the initial stress remains.

For a one-dimensional Maxwell model:
$$
\tau_M=\frac{\eta}{E}
$$
$\tau_M$ = relaxation time (s); $\eta$ = dashpot viscosity (Pa s); $E$ = spring's Young's modulus (Pa). A dashpot represents viscous resistance to deformation. Real materials can have several relaxation times; this simple model is not universal.

**Homologous temperature**
Temperature relative to a material's melting temperature; helps indicate how readily time-dependent deformation and structural recovery can occur.

$$
T_H=\frac{T}{T_m}
$$
$T_H$ = homologous temperature (dimensionless); $T_m$ = melting temperature. Both temperatures must be in kelvin. For lithium at 25°C, $T_H\approx298/454\approx0.66$.

**Annealing**
Thermally activated changes in a material's internal structure that reduce stored energy, such as defects rearranging or disappearing. It usually involves heating, but lithium can undergo recovery even at room temperature because its homologous temperature is high. The material remains solid.

![[Attachments/Concept diagrams/annealing-grain-structure.png|700]]

*Rolling elongates grains and introduces defects. Annealing can reduce defects, form new grains (recrystallisation), and allow grains to grow. This is a general metal-processing illustration; lithium's room-temperature recovery need not include all these stages.*

**Recrystallisation** — Formation of new grains with fewer deformation-induced defects.
**Grain growth** — Larger grains grow at the expense of smaller grains.

## Mechanics: modelling

**Homogeneous**
Uniform in composition or specified properties from place to place at the scale considered. Different from isotropic, which means the same in every direction.

**Heterogeneous**
Having composition or properties that vary from place to place, such as a composite cathode containing different materials.

**Heterogeneity**
Spatial variation in composition, structure, or properties; for example, uneven porosity or stiffness.

**Heterogenisation**
The process of making a material or field less uniform, so differences develop from place to place.

**Homogenisation**
A modelling method that replaces detailed small-scale variations with effective, averaged properties at a larger scale.

**Continuum**
A description that treats material as continuously distributed, with fields such as stress and concentration, rather than tracking individual atoms. It can still include spatial variations and interfaces.

**Constitutive law**
An equation describing material behaviour, such as the relationship between stress and strain.

For small-strain linear elasticity with chemical expansion:
$$
\boldsymbol{\sigma}=\mathbb{C}:\boldsymbol{\varepsilon}^{e}
$$
$\mathbb{C}$ = elastic stiffness tensor; $\boldsymbol{\sigma}$ = stress tensor; $\boldsymbol{\varepsilon}^{e}$ = elastic strain tensor. The colon means multiplying the stiffness components by the strain components and summing them.

**Mechanical equilibrium**
A state in which forces balance.

$$
\nabla\cdot\boldsymbol{\sigma}+\mathbf{b}=\mathbf{0}
$$
$\boldsymbol{\sigma}$ = stress tensor; $\mathbf{b}$ = body force per unit volume, such as gravity (N/m$^3$); $\nabla\cdot$ = divergence, measuring the net internal force from spatial stress variation. Inertia is neglected.

**Quasistatic model**
A model that assumes forces balance as deformation evolves and neglects inertia.

**Traction**
Force per unit area acting on a surface.

$$
\mathbf{t}=\boldsymbol{\sigma}\mathbf{n}
$$
$\mathbf{t}$ = traction vector (force per unit area); $\boldsymbol{\sigma}$ = stress tensor; $\mathbf{n}$ = outward unit normal to the surface. A free surface has $\mathbf{t}=\mathbf{0}$.

**Mechanical boundary conditions**
The imposed movements or forces at a model's boundaries.

**Small deformation**
Deformation with small strains and rotations, allowing simplified equations.

$$
\boldsymbol{\varepsilon}=\frac{1}{2}\left(\nabla\mathbf{u}+(\nabla\mathbf{u})^T\right)
$$
$\boldsymbol{\varepsilon}$ = small-strain tensor; $\mathbf{u}$ = displacement vector; $\nabla\mathbf{u}$ = displacement gradient; $T$ = matrix transpose (swap rows and columns). This approximation requires small displacement gradients.

**Finite deformation**
A description that accounts for changes in geometry, including large strains or rotations.

## Mechanics: interfaces and damage

**Contact mechanics**
The study of forces and deformation where surfaces touch.

**Contact loss**
Separation of surfaces that previously touched, which can interrupt transport or reactions.

**Stack pressure**
External compression applied across a battery cell.

$$
p_{\mathrm{stack}}=\frac{F_{\mathrm{compression}}}{A_{\mathrm{cell}}}
$$
$F_{\mathrm{compression}}$ = applied compressive force; $A_{\mathrm{cell}}$ = loaded cell area. This is an overall applied pressure, not the stress at each particle contact.

**Decohesion / delamination**
Separation along an interface between materials or layers.

**Stress concentration**
A local increase in stress near a crack, pore, or sharp feature.

**Strength**
The stress a material can withstand before a specified failure, such as yielding or breaking.

**Fracture**
Cracking or breaking of a material.

**Fracture toughness**
A material's resistance to the growth of a crack.

**Pugh's ratio**
A ratio of elastic stiffnesses used as a rough indicator of ductile versus brittle behaviour.

Using the bulk-to-shear convention:
$$
\text{Pugh's ratio}=\frac{B}{G}
$$
$B$ = bulk modulus: resistance to elastic volume change under uniform pressure; $G$ = shear modulus: resistance to elastic shape change under shear. Both have units of Pa, so the ratio is dimensionless.

A common rule of thumb is $B/G\gtrsim1.75$ suggests more ductile behaviour, while lower values suggest more brittle behaviour. This is an empirical guide, not a universal boundary or a measurement of fracture toughness. Some papers use the inverse $G/B$, so always check the definition.

Source: [Example research paper explaining the criterion](https://pmc.ncbi.nlm.nih.gov/articles/PMC5571045/).

**Cohesive energy**
The energy needed to separate a solid into isolated atoms; a measure of how strongly it is bound together. Higher cohesive energy does not necessarily mean higher fracture toughness: a strongly bonded material can still be brittle if it cannot deform plastically near a crack.

**Fatigue**
Damage caused by repeated loading and unloading.

## Crystal structure

### Basic concepts

**Crystal**
A solid with an ordered atomic structure; ordinary crystals have a repeating pattern.

**Crystal lattice**
A repeating grid of equivalent points describing a crystal's periodicity.

**Basis (atomic motif)**
The atom or group of atoms repeated at each lattice point.

**Crystal structure**
The arrangement of atoms formed by combining a lattice with its atomic basis.

**Unit cell**
A repeating building block that reproduces the crystal when translated in three directions.

**Primitive unit cell**
A smallest-volume unit cell containing one lattice point, counting shared points fractionally.

**Conventional unit cell**
A unit cell chosen to make the crystal's symmetry easy to see.

**Lattice parameter**
A length or angle defining the repeating unit of a crystal's atomic structure.

Relative change along one lattice direction:
$$
\varepsilon_a=\frac{a-a_0}{a_0}
$$
$a_0$ = reference lattice length; $a$ = current lattice length; $\varepsilon_a$ = relative lattice change. It can include both chemical and elastic contributions.

**Lattice spacing**
The distance between neighbouring parallel planes in a crystal.

**Unit-cell volume**
The volume of the crystal's repeating unit, determined by its lattice parameters.

For a cubic cell:
$$
V_{\mathrm{cell}}=a^3
$$
$a$ = cube edge length; $V_{\mathrm{cell}}$ = unit-cell volume. For small equal changes in all three directions, $\Delta V_{\mathrm{cell}}/V_{\mathrm{cell},0}\approx3\Delta a/a_0$, where $\Delta$ means change and subscript $0$ means the reference state.

**Crystal system**
A classification by crystal symmetry, such as cubic, hexagonal, or tetragonal.

**Bravais lattice**
One of the 14 types of repeating three-dimensional lattices.

### Common structures
![[Pasted image 20261009114953.png]]

**Simple cubic (SC)**
A cubic lattice with points at the cube's corners only.

**Body-centred cubic (BCC)**
A cubic lattice with points at the corners and one at the cube's centre.

**Face-centred cubic (FCC)**
A cubic lattice with points at the corners and the centre of each face.

**Hexagonal close-packed (HCP)**
A structure of closely packed atomic layers stacked in an ABAB pattern.

**Close packing**
An arrangement that packs equal-sized spheres as densely as possible; FCC and HCP are examples.

**Coordination number**
The number of nearest neighbours around an atom or ion.

**Layered structure**
A crystal structure with distinct atomic layers, often allowing lithium to move between them.

**Garnets (garnet-type materials)**
A family of materials with a related crystal structure. In SSBs, lithium-conducting garnet-type ceramics such as LLZO ($\mathrm{Li_7La_3Zr_2O_{12}}$) are used as solid electrolytes: lithium ions move between sites through the solid crystal. Not all garnets are good ion conductors.

Source: [LLZO structure and ion transport study](https://www.nature.com/articles/srep18053).

### Grains, directions, and defects

![[Attachments/Concept diagrams/single-polycrystalline-amorphous.png|700]]

*Left: single crystal* — one ordered crystal orientation. *Middle: polycrystalline* — ordered grains with different orientations, separated by grain boundaries. *Right: amorphous* — no long-range repeating atomic pattern, although local order can exist. The left label “Crystalline” means single-crystalline in this sketch; both left and middle are crystalline. Circles represent atoms schematically.

Image supplied by you; original source not identified.

**Single crystal**
A crystal with one continuous crystallographic orientation throughout.

**Polycrystalline material**
A material made of many crystal grains, often with different orientations.

**Grain**
A region of a material with a continuous crystal structure and orientation.

**Grain boundary**
The interface between crystal grains with different orientations.

![[Attachments/Concept diagrams/grain-boundary.jpg|314]]

*Grains and boundaries:* each numbered region is a grain; the lines between them are grain boundaries. Here $d$ means a characteristic grain size. This sketch shows grain regions, not individual atoms or battery particles.

Image: Wlytle1375, [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Grain_boundary_diagram.jpg), public domain. Image unchanged.

**Amorphous material**
A material without long-range periodic atomic order.

**Miller indices**
Numbers such as (100) or (111) identifying the orientation of crystal planes.

**Vacancy**
A missing atom or ion at a normally occupied crystal site.

**Interstitial site**
A position between the usual atomic sites where an atom or ion can sit.

**Dislocation**
A line defect in a crystal whose movement can produce plastic deformation.
![[Pasted image 20261009164747.png]]

**Shear-driven dislocation glide**
Movement of a dislocation along a crystal plane under shear stress, producing permanent slip (plastic deformation). Like moving a wrinkle across a rug, atoms rearrange locally rather than a whole plane sliding at once. This can help relieve stress when dislocations can move.

**Dislocation nucleation–governed plasticity**
Plastic deformation controlled mainly by creating new dislocations, which then glide and produce permanent slip. When few mobile dislocations are available, forming new ones can require high stress. This is one possible explanation for the high strength of very small lithium volumes.

*Nucleation* means formation of a new defect. Here it refers to dislocations, not the formation of new lithium deposits.

Source: [Small-scale lithium strength study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5224391/).

**Slip plane**
The crystal plane along which a dislocation glides.

**Burgers vector**
A vector describing the size and direction of the crystal slip associated with a dislocation.

![[Attachments/Concept diagrams/edge-dislocation.jpg|350]]

*An edge dislocation:* an extra half-plane ends inside the crystal; its terminating edge is the blue dislocation line. Under suitable shear stress, this defect can glide through the crystal and leave permanent slip. $b$ labels the Burgers vector. This image shows the defect, rather than a sequence of its motion.

For a straight dislocation, the magnitude of the glide force per unit length is:
$$
f_g=|\tau|b
$$
$f_g$ = glide force per unit dislocation length (N/m); $\tau$ = shear stress on the slip plane in the slip direction (Pa); $b$ = Burgers-vector magnitude (m). Greater shear stress gives a greater driving force; obstacles can still resist motion.

Explanation and equation: [Princeton materials teaching notes](https://www.princeton.edu/~maelabs/mae324/07/07mae_23.htm), [dislocations and slip](https://www.princeton.edu/~maelabs/hpt/materials/dislocations.htm).
Image: Wikityke, [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Dislocation_edge_d2.jpg), [CC BY-SA 2.5](https://creativecommons.org/licenses/by-sa/2.5/). Image unchanged.

## Coupling mechanics and chemistry

**Chemical strain (eigenstrain)**
The expansion or contraction a material would undergo freely when its composition changes.

For small strain, ignoring thermal expansion:
$$
\boldsymbol{\varepsilon}=\boldsymbol{\varepsilon}^{e}+\boldsymbol{\varepsilon}^{\mathrm{chem}}+\boldsymbol{\varepsilon}^{p}
$$
The terms are total, elastic, chemical, and plastic strain tensors. Stress depends on the elastic part, not directly on the total strain.

A simple isotropic, linear expansion model is:
$$
\boldsymbol{\varepsilon}^{\mathrm{chem}}=\frac{\Omega(c-c_0)}{3}\mathbf{I}
$$
$c$ = lithium concentration (mol/m$^3$); $c_0$ = reference concentration; $\Omega$ = constant chemical expansion coefficient (m$^3$/mol), often modelled as lithium's partial molar volume; $\mathbf{I}$ = identity tensor. Division by three shares the volume change equally among three directions. Layered materials may need an anisotropic model instead.

**Reversible strain**
Deformation that disappears when its cause is reversed. This includes elastic strain and reversible chemical strain, such as expansion during lithium insertion that reverses during removal.

**Chemical potential**
The free-energy change per mole of lithium added; differences in chemical potential drive diffusion.

A simple small-strain model for lithium in an isotropically expanding host is:
$$
\mu=\mu_{\mathrm{chem}}(c)-\Omega\sigma_h
$$
$\mu$ = lithium chemical potential (J/mol); $\mu_{\mathrm{chem}}(c)$ = chemical contribution at concentration $c$; $\Omega$ = constant chemical expansion coefficient (m$^3$/mol); $\sigma_h$ = hydrostatic stress, positive in tension. This simplified expression omits other contributions and is not a general ion-transport law for the electrolyte.

**Diffusion length (for diffusional stress relief)**
The characteristic distance atoms must travel to redistribute material and accommodate deformation. The size and shape of the stressed lithium region help determine this distance; it is not necessarily the whole electrolyte thickness.

A simple diffusion-time estimate is:
$$
t_{\mathrm{diff}}\sim\frac{L^2}{D}
$$
$t_{\mathrm{diff}}$ = diffusion time (s); $L$ = characteristic travel distance (m); $D$ = relevant diffusion coefficient (m$^2$/s); $\sim$ means an order-of-magnitude estimate. Doubling $L$ roughly quadruples this time. This estimates redistribution time, not a complete creep law; geometry and diffusion route also matter.

**Stressed volume**
The region of material experiencing stress. Its dimensions and shape influence diffusion distances and the number of available dislocations.

**Diffusion flux**
The amount of lithium crossing a unit area per unit time, with a direction.

For an isotropic mobility model in an active material:
$$
\mathbf{J}=-\frac{Dc}{RT}\nabla\mu
$$
$\mathbf{J}$ = lithium molar flux (mol/m$^2$/s); $D$ = mobility-related diffusivity (m$^2$/s); $c$ = lithium concentration; $R$ = gas constant; $T$ = absolute temperature; $\nabla\mu$ = spatial gradient of chemical potential. The minus sign gives flow towards lower chemical potential. A stress gradient can therefore affect diffusion. This $D$ need not equal a measured chemical diffusion coefficient.

## Electrochemistry

**Redox reaction**
A reaction involving electron transfer: one species loses electrons (oxidation), while another gains them (reduction).

Lithium deposition at a lithium-metal electrode:
$$
\mathrm{Li^+}+e^-\rightleftharpoons\mathrm{Li\,(metal)}
$$
$e^-$ = electron; the rightward reaction deposits lithium, while the reverse removes it. This equation describes lithium metal, not insertion into a cathode host.
