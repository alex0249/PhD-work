
# Workflow

 Author Year — Short title

Reference / DOI:
PDF link:
Status: skimmed / studying / read

**Why I am reading this**
Which question in my PhD does it help answer?

**Main contribution**
Explain in 3–5 sentences, in my own words.

**Model or main details**
- Geometry and materials:
- Unknown variables:
- Governing equations:
- Initial and boundary conditions:
- Main assumptions:
- Coupling between mechanics and electrochemistry:

**Equations I need**
Equation number and page.
Meaning of each term, units, and when it applies.

**Evidence**
What do they compare with: experiments, analytical
solutions, or another numerical model?
What does that comparison actually establish?

**Relevance and limitations**
What could I use?
Which assumptions might fail for my battery?

**Questions and next action**
What do I not understand?
What will I derive, check, implement, or ask?

# Kalnaus review of mechanics in SSBs
[Link](file:///C:/Users/alexb/AppData/Local/Temp/d6cc8fd7-3f96-46db-9b90-758f3f308f19_OneDrive_1_10-9-2026.zip.f19/2023%20Kalnaus%20Solid-state%20batteries%20the%20critical%20role%20of%20mechanics.pdf)
**Why:** 
Overview of why mechanics is important in SSBs

**Main details:** 
- Clear summary of mechanical effects that need solving: (i) plating on solid electrolyte surface (ii) contact loss at electrode contacts and grain boundaries from volume changes during cycling (iii) manufacturing difficulties
- Lithium metal deposition on surface + volume defects of ceramic electrolyte  => local high stresses => electrolyte fracture + metallic lithium in cracks
- Manufacturing: cathode-electrolyte stacks need to be able to withstand forces from equipment

**Details**
Cathode (layered oxide): 
- Li+ creates gradient in lattice parameter = distances change with position inside the particle as: lithium changes crystal dimensions when diffusing through particle
- Deviatoric strain + dilatation: distortion of particle shape and volume change
	=> lithium entering particle does not necessarily expand equally in every direction
	=> if tensile component of interfacial stress > strength of stressed material => brittle fracture
	
Anode (Li metal): 
- Li depos non-uniform => localised stresses at interface => unstable and roughen Li => dendrites/filaments
- Sustainable pressure in small, constrained volumes of Li is dep. on strain rate which dep. on oper. current dens.  

Needed mechanical capabilities during rapid volume changes:
- minimise localised gradients in elastic strain
- efficient stress relief mechanisms that work over many operating conditions
- no loss of contact => may need a form of reversible strain

Zero-strain (ZS) cathodes:
- no change in lattice params during cycling
- current materials are too low V 
- NMC can be dev. that min. vol. changes as fn. of SoC
- chem comp, ionic ordering & metal coord. are critical for ZS cath

Pairing with neg strain mat with pos strain mat:
- balance vol exp and contr 
- mitigate mech failures at cath-electrolyte interface but local strain + stress may still lead to deg

Planar plating and stripping of Li metal
- optimal to enable planar layer-by-layer, stress-free plating and stripping of Li metal
	=> difficult as non-uni dep is not well understood for single-ion conducting solid electrolytes
- improving homogeneity of Li plating and stripping is v beneficial: Ag-C interface layer, improving wetting
- cond of Li is important => defects = bad plating

Stress-relief mechanisms:
- Deformation can relieve stress: goal is to activate inelastic or viscoelastic strain to reduce stress mag
- Diff mat relieve stress differently
- Mech must act fast enough: high I => faster def => less time for stress relief
- Temp and length scale matter: mech that works at one condition might not work at another
- If stress cannot relax sufficiently => cracking can occur
Question: Why is plastic deformation a good stress relief mechanism? 

Plastic deformation of Li metal
- bulk lithium is soft, but lithium in very small, constrained regions can withstand much higher stresses.
- Li is close to melting temp at room temp => promotes thermally activated glide
- Li metal continuously annealing 
=> net effect: bulk Li has no notable capac for work hardening at room temp unless strain rate is suff high
	=>Li can accom shape ch and contact with electrolyte as Li is removed
	=> easy flow => push Li into electrolyte pores or defects => short circuit
- Small, constrained volumes behave differently: exp show that plastic def in bulk Li metal limits flow stress to 1MPa or less but within constrained volumes of Li metal oper effic of dislocation glide can become severely limited even at high T
- Li trapped in electrolyte defects develop much higher stresses than bulk 
- Diffusion is proposed as the dom mech of action for deformation
- "Danger zone" proposed which can form at critical length scale
- Narayan used large deformation theory to model shear flow of lithium: suggested the incorporation of Nabarro-Herring (volume diffusion), HarperDorn (nonconservative dislocation motion), and interface diffusional flow 
- also need to understand the effect of external stack pressure 

Fracture: 
- stress can be transferred to the SSE: if Li cannot deform sufficiently or cathode is non-uniform which leads to localised gradients in elastic strain
- ceramic electrolytes and cathodes are usually brittle: dislocation nucleation requires large energy penalty => material cant move/change shape so it cracks 
- Fracture toughness, $K_C$: measures resistance to fracture- Low toughness means cracks can grow relatively easily; the required stress also depends on crack size and geometry. 
	=> ionic crystals (cathodes and ceramic electrolyte) are not that tough => will crack
- removing grain boundaries in LLZO did not change resistance to cracking or its low resistance (LLZO is a ceramic garnet SSE material)
- Need materials that produce no hard contacts between SSE and cathode (concentrated stress) + can relieve stress without fracturing (allow for plastic deformation)

Plastic deformation in ceramics:
- three sources of stress in SSB: (i) Li plating into defects in solid electrolyte (ii) constrained cathode expansion (iii) external pressure
- ceramics don't deform plastically at room temperature => crack
- one approach: introduce defects/dislocations during manufacturing => move existing dislocations => if high dens of disl near crack tip => plastic deform => reduce stress conc
		- has been shown experimentally
- noncrystalline electrolytes (e.g. glasses) have different deform mech: despite being brittle at macroscale,  local densification and shear flow can relieve stress without disl motion => LiPON is promising

Time-dependent mechanical behaviour:
- stress-induced densification in glass can be partially recoverable unlike permanent plastic shear flow => annealing at temp < glass trans temp
- silica is an anomalous glass as inelastic deformation occurs just through densification
- Glasses usually flow extremely slowly below their glass-transition temperature, but high local stresses can activate measurable creep and viscoelastic deformation even at room temperature.
- Repeated loading and unloading can produce hysteresis: the material follows different deformation paths during loading and recovery.
