
 .. role:: cpp(code)
    :language: c++

.. _Equations:



Equations
=========

TODO
----

placeholder

.. math::
 
  \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u}) + S_{{\rm ext},\rho},

  \frac{\partial (\rho \mathbf{u})}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u} + \mathbf{\Pi}) - \nabla p +\rho \mathbf{g} + \mathbf{S}_{{\rm ext},\rho\mathbf{u}},

  \frac{\partial (\rho E)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} E + p \mathbf{u}) - \nabla \cdot (\mathbf{\Pi} \cdot \mathbf{u}) + \rho \mathbf{u} \cdot \mathbf{g} + \nabla\cdot \boldsymbol{\mathcal{Q}}+ S_{{\rm ext},\rho E},

  \frac{\partial (\rho Y_k)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} Y_k)
  - \nabla \cdot \boldsymbol{\mathcal{F}}_{k} + \rho \dot\omega_k + S_{{\rm ext},\rho Y_k},

  \frac{\partial (\rho A_k)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} A_k) + S_{{\rm ext},\rho A_k},

  \frac{\partial (\rho B_k)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} B_k) + S_{{\rm ext},\rho B_k}.


Here :math:`\rho, \mathbf{u}, T`, and :math:`p` are the density, velocity,
temperature and pressure, respectively. :math:`E
= e + \mathbf{u} \cdot \mathbf{u} / 2` is the total energy with :math:`e` representing the internal energy, which is defined as in the CHEMKIN standard to include both sensible
and chemical energy (species heats of formation) and is conserved across chemical reactions. 
:math:`Y_k` is the mass fraction of the :math:`k^{\rm th}` species,
with associated production rate, :math:`\dot\omega_k`.  Here :math:`\mathbf{g}` is the gravitational vector, and
:math:`S_{{\rm ext},\rho}, \mathbf{S}_{{\rm ext},\rho\mathbf{u}}`, etc., are user-specified
source terms.  :math:`A_k` is an advected quantity, i.e., a tracer.  Also
:math:`\boldsymbol{\mathcal{F}}_{m}, \mathbf{\Pi}`, and :math:`\boldsymbol{\mathcal{Q}}` are
the diffusive transport fluxes for species, momentum and heat.  Note that the internal
energy for species :math:`k` includes its heat of formation (and can therefore take on negative and
positive values). The auxiliary fields, :math:`B_k`, have a user-defined
evolution equation, but by default are treated as advected quantities.
