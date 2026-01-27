"""
Unit tests for axisymmetric physics implementation.
"""
import pytest
import numpy as np
from dolfin import UnitSquareMesh, Function, VectorFunctionSpace, FunctionSpace, Expression

# Import physics functions
from src.core.physics import (
    epsilon, epsilon_axi,
    sigma, sigma_axi,
    psi_positive, psi_positive_axi
)


class TestAxiStrain:
    """Tests for axisymmetric strain tensor."""
    
    def test_epsilon_axi_shape(self):
        """Verify epsilon_axi returns 3x3 tensor."""
        mesh = UnitSquareMesh(4, 4)
        V = VectorFunctionSpace(mesh, "CG", 1)
        u = Function(V)
        
        from dolfin import SpatialCoordinate
        x = SpatialCoordinate(mesh)
        r = x[0]
        
        eps = epsilon_axi(u, r)
        
        # Check it's a 3x3 tensor (UFL shape)
        assert eps.ufl_shape == (3, 3), f"Expected (3,3), got {eps.ufl_shape}"
    
    def test_epsilon_2d_shape(self):
        """Verify standard epsilon returns 2x2 tensor."""
        mesh = UnitSquareMesh(4, 4)
        V = VectorFunctionSpace(mesh, "CG", 1)
        u = Function(V)
        
        eps = epsilon(u)
        assert eps.ufl_shape == (2, 2), f"Expected (2,2), got {eps.ufl_shape}"


class TestAxiStress:
    """Tests for axisymmetric stress tensor."""
    
    def test_sigma_axi_shape(self):
        """Verify sigma_axi returns 3x3 tensor."""
        mesh = UnitSquareMesh(4, 4)
        V = VectorFunctionSpace(mesh, "CG", 1)
        u = Function(V)
        
        from dolfin import SpatialCoordinate
        x = SpatialCoordinate(mesh)
        r = x[0]
        
        lmbda = 1e8
        mu = 5e7
        
        sig = sigma_axi(u, r, lmbda, mu)
        assert sig.ufl_shape == (3, 3), f"Expected (3,3), got {sig.ufl_shape}"
    
    def test_sigma_2d_shape(self):
        """Verify standard sigma returns 2x2 tensor."""
        mesh = UnitSquareMesh(4, 4)
        V = VectorFunctionSpace(mesh, "CG", 1)
        u = Function(V)
        
        lmbda = 1e8
        mu = 5e7
        
        sig = sigma(u, lmbda, mu)
        assert sig.ufl_shape == (2, 2), f"Expected (2,2), got {sig.ufl_shape}"


class TestEnergyDensity:
    """Tests for strain energy density functions."""
    
    def test_psi_positive_returns_scalar(self):
        """Verify psi_positive returns scalar form."""
        mesh = UnitSquareMesh(4, 4)
        V = VectorFunctionSpace(mesh, "CG", 1)
        u = Function(V)
        
        lmbda = 1e8
        mu = 5e7
        
        psi = psi_positive(u, lmbda, mu)
        # Scalar has empty shape
        assert psi.ufl_shape == (), f"Expected (), got {psi.ufl_shape}"
    
    def test_psi_positive_axi_returns_scalar(self):
        """Verify psi_positive_axi returns scalar form."""
        mesh = UnitSquareMesh(4, 4)
        V = VectorFunctionSpace(mesh, "CG", 1)
        u = Function(V)
        
        from dolfin import SpatialCoordinate
        x = SpatialCoordinate(mesh)
        r = x[0]
        
        lmbda = 1e8
        mu = 5e7
        
        psi = psi_positive_axi(u, r, lmbda, mu)
        assert psi.ufl_shape == (), f"Expected (), got {psi.ufl_shape}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
