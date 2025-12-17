# Topological Shape Detection in Julia

## Overview

This project explores **topological data analysis (TDA)** techniques to **detect, distinguish, and classify geometric shapes** from sampled data such as point clouds.  
The implementation is entirely done in **Julia**, taking advantage of its performance and numerical capabilities.

The main idea is that **topology captures global geometric features** (connected components, holes, voids) that are robust to noise, rotations, and translations, making it well-suited for shape recognition.

---

## Motivation

In many real-world datasets (biology, medical imaging, physics, computer vision), shapes are not explicitly given but must be inferred from data.  
Classical geometric methods often fail when shapes are noisy, deformed, or embedded in higher dimensions.

**Topological methods**, such as homology and persistent homology, provide invariants that:
- are stable under perturbations,
- do not depend on exact coordinates,
- can distinguish shapes that look similar geometrically but differ topologically.

---

## Shapes Studied

We focus on recognizing and distinguishing shapes such as:
- Line segment  
- Circle  
- Disk  
- Sphere  
- Torus  
- Ellipsoid  
- Irregular (potato-like) 3D shapes  

Each shape is represented by a **point cloud** or a **cellular complex**.

---

## Methodology

The project follows these main steps:

1. **Shape generation / loading**  
   - Generate synthetic point clouds for known shapes  
   - Load provided datasets in higher dimensions (e.g. ℝ⁴)

2. **Topological approximation**  
   - Build simplicial or cubical complexes  
   - Use filtrations (e.g. Vietoris–Rips complexes)

3. **Topological feature extraction**  
   - Compute homology groups  
   - Compute **persistent homology**  
   - Extract persistence diagrams and derived representations

4. **Shape discrimination**  
   - Compare persistence diagrams  
   - Use topological signatures (e.g. Betti numbers, persistence silhouettes)

5. **Robustness testing**  
   - Apply rotations and translations  
   - Check invariance of topological descriptors

---

## Example

### Distinguishing a Circle from a Disk

- A **circle** has:
  - 1 connected component  
  - 1 persistent 1-dimensional hole  

- A **disk** has:
  - 1 connected component  
  - no persistent 1-dimensional hole  

Using persistent homology:
- The circle produces a **long-lived 1D feature** in the persistence diagram
- The disk does not

This allows us to reliably distinguish the two shapes even with noisy samples.

---

## Technologies

- **Language**: Julia  
- **Paradigm**: Scientific computing, computational topology  
- **Key concepts**:
  - Simplicial complexes
  - Persistent homology
  - Topological invariants



## Activate the environnement 

 using Pkg
Pkg.activate(".")

To verify you are on the environment use:
Pkg.status()


## run the test 
import Pkg
Pkg.activate(".")
Pkg.test()