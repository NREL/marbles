#ifndef STENCIL_H
#define STENCIL_H

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include "Constants.H"

namespace lbm::stencil {
struct Stencil
{
    Stencil() = default;

#if AMREX_SPACEDIM == 2
    //Warning! Achtung! This stencil is for 3D only. 
    //Make Z direction periodic and keep the number of points 1 in Z direction for 2D simulation
    //For a strictly 2D lattice, please use the Stencil.H header

    const amrex::GpuArray<amrex::IntVect, constants::N_MICRO_STATES> evs = {
        amrex::IntVect(0, 0),   amrex::IntVect(1, 0),  amrex::IntVect(-1, 0),
        amrex::IntVect(0, 1),   amrex::IntVect(0, -1), amrex::IntVect(1, 1),
        amrex::IntVect(-1, -1), amrex::IntVect(1, -1), amrex::IntVect(-1, 1)};

    const amrex::GpuArray<int, constants::N_MICRO_STATES> bounce_dirs = {
        0, 2, 1, 4, 3, 6, 5, 8, 7};

    const amrex::GpuArray<amrex::Real, constants::N_MICRO_STATES> weights = {
        4.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0, 1.0 / 9.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};

    //const amrex::GpuArray<int, constants::N_MICRO_STATES> bounce_Xdir={};
    //const amrex::GpuArray<int, constants::N_MICRO_STATES> bounce_Ydir={};
    //const amrex::GpuArray<int, constants::N_MICRO_STATES> bounce_Zdir={};
    

#else
    // const amrex::GpuArray<amrex::IntVect, constants::N_MICRO_STATES> evs = {
    //     amrex::IntVect(0, 0, 0),    amrex::IntVect(1, 0, 0),
    //     amrex::IntVect(-1, 0, 0),   amrex::IntVect(0, 1, 0),
    //     amrex::IntVect(0, -1, 0),   amrex::IntVect(0, 0, 1),
    //     amrex::IntVect(0, 0, -1),   amrex::IntVect(1, 1, 0),
    //     amrex::IntVect(-1, -1, 0),  amrex::IntVect(1, -1, 0),
    //     amrex::IntVect(-1, 1, 0),   amrex::IntVect(1, 0, 1),
    //     amrex::IntVect(-1, 0, -1),  amrex::IntVect(1, 0, -1),
    //     amrex::IntVect(-1, 0, 1),   amrex::IntVect(0, 1, 1),
    //     amrex::IntVect(0, -1, -1),  amrex::IntVect(0, 1, -1),
    //     amrex::IntVect(0, -1, 1),   amrex::IntVect(1, 1, 1),
    //     amrex::IntVect(-1, -1, -1), amrex::IntVect(1, 1, -1),
    //     amrex::IntVect(-1, -1, 1),  amrex::IntVect(1, -1, 1),
    //     amrex::IntVect(-1, 1, -1),  amrex::IntVect(-1, 1, 1),
    //     amrex::IntVect(1, -1, -1)};

    // const amrex::GpuArray<int, constants::N_MICRO_STATES> bounce_dirs = {
    //     0,  2,  1,  4,  3,  6,  5,  8,  7,  10, 9,  12, 11, 14,
    //     13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};

    // const amrex::GpuArray<amrex::Real, constants::N_MICRO_STATES> weights = {
    //     8.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,
    //     2.0 / 27.0,  2.0 / 27.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,
    //     1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,
    //     1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 216.0,
    //     1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0,
    //     1.0 / 216.0, 1.0 / 216.0};

    const int cSC1 = 1; //simple cubic
    const int cFCC1= 1; //face centered cubic
    const int cBCC1= 1; //body centered cubic

    const amrex::GpuArray<amrex::IntVect, constants::N_MICRO_STATES> evs = {
        amrex::IntVect( 0, 0, 0),

        amrex::IntVect( cSC1, 0, 0), amrex::IntVect(-cSC1, 0, 0),
        amrex::IntVect( 0, cSC1, 0), amrex::IntVect( 0,-cSC1, 0),
        amrex::IntVect( 0, 0, cSC1), amrex::IntVect( 0, 0,-cSC1),

        amrex::IntVect( cFCC1, cFCC1, 0), amrex::IntVect(-cFCC1,-cFCC1, 0),  
        amrex::IntVect( cFCC1,-cFCC1, 0), amrex::IntVect(-cFCC1, cFCC1, 0),
        amrex::IntVect( cFCC1, 0, cFCC1), amrex::IntVect(-cFCC1, 0,-cFCC1),
        amrex::IntVect( cFCC1, 0,-cFCC1), amrex::IntVect(-cFCC1, 0, cFCC1),
        amrex::IntVect( 0, cFCC1, cFCC1), amrex::IntVect( 0,-cFCC1,-cFCC1),
        amrex::IntVect( 0, cFCC1,-cFCC1), amrex::IntVect(0, -cFCC1, cFCC1),

        amrex::IntVect( cBCC1, cBCC1, cBCC1), amrex::IntVect(-cBCC1,-cBCC1,-cBCC1),
        amrex::IntVect( cBCC1,-cBCC1, cBCC1), amrex::IntVect(-cBCC1, cBCC1,-cBCC1),
        amrex::IntVect( cBCC1,-cBCC1,-cBCC1), amrex::IntVect(-cBCC1, cBCC1, cBCC1),        
        amrex::IntVect( cBCC1, cBCC1,-cBCC1), amrex::IntVect(-cBCC1,-cBCC1, cBCC1)  
        };

        const amrex::Real wSC1   = 2.0/27.0;
        const amrex::Real wFCC1  = 1.0/54.0;
        const amrex::Real wBCC1  = 1.0/216.0;
        const amrex::Real w0 = 1.0 - (6.0*(wSC1) + 12.0*(wFCC1) + 8.0*(wBCC1));

        const amrex::GpuArray<amrex::Real, constants::N_MICRO_STATES> weights = {
        	w0,	
            wSC1,	wSC1,	wSC1,	wSC1,	wSC1,	wSC1,
            wFCC1,	wFCC1,	wFCC1,	wFCC1,	wFCC1,	wFCC1,	wFCC1,	wFCC1,	wFCC1,	wFCC1,	wFCC1,	wFCC1,
            wBCC1,	wBCC1,	wBCC1,	wBCC1,	wBCC1,	wBCC1,	wBCC1,	wBCC1 
            };

        enum dvD3Q27{
            DV_ZERO_ZERO_ZERO, 
            
            DV_P1_ZERO_ZERO,DV_M1_ZERO_ZERO,
            DV_ZERO_P1_ZERO,DV_ZERO_M1_ZERO,
            DV_ZERO_ZERO_P1,DV_ZERO_ZERO_M1,
            
            DV_P1_P1_ZERO,DV_M1_M1_ZERO,
            DV_P1_M1_ZERO,DV_M1_P1_ZERO,
            DV_P1_ZERO_P1,DV_M1_ZERO_M1,
            DV_P1_ZERO_M1,DV_M1_ZERO_P1,
            DV_ZERO_P1_P1,DV_ZERO_M1_M1,
            DV_ZERO_P1_M1,DV_ZERO_M1_P1,

            DV_P1_P1_P1,DV_M1_M1_M1,
            DV_P1_M1_P1,DV_M1_P1_M1,
            DV_P1_M1_M1,DV_M1_P1_P1,
            DV_P1_P1_M1,DV_M1_M1_P1
            };

        const amrex::GpuArray<int, constants::N_MICRO_STATES> bounce_dirs = {
        DV_ZERO_ZERO_ZERO,
        DV_M1_ZERO_ZERO,DV_P1_ZERO_ZERO,
        DV_ZERO_M1_ZERO,DV_ZERO_P1_ZERO,
        DV_ZERO_ZERO_M1,DV_ZERO_ZERO_P1,

        DV_M1_M1_ZERO,DV_P1_P1_ZERO,
        DV_M1_P1_ZERO,DV_P1_M1_ZERO,
        DV_M1_ZERO_M1,DV_P1_ZERO_P1,
        DV_M1_ZERO_P1,DV_P1_ZERO_M1,
        DV_ZERO_M1_M1,DV_ZERO_P1_P1,
        DV_ZERO_M1_P1,DV_ZERO_P1_M1,

        DV_M1_M1_M1,DV_P1_P1_P1,
        DV_M1_P1_M1,DV_P1_M1_P1,
        DV_M1_P1_P1,DV_P1_M1_M1,
        DV_M1_M1_P1,DV_P1_P1_M1
        };

    const amrex::GpuArray<int, constants::N_MICRO_STATES> bounce_Xdir={
    DV_ZERO_ZERO_ZERO, 
    DV_M1_ZERO_ZERO,   DV_P1_ZERO_ZERO,   
    DV_ZERO_P1_ZERO,   DV_ZERO_M1_ZERO,   
    DV_ZERO_ZERO_P1,   DV_ZERO_ZERO_M1,   

    DV_M1_P1_ZERO,   DV_P1_M1_ZERO,    
    DV_M1_M1_ZERO,   DV_P1_P1_ZERO,   
    DV_M1_ZERO_P1,   DV_P1_ZERO_M1,   
    DV_M1_ZERO_M1,   DV_P1_ZERO_P1,   
    DV_ZERO_P1_P1,   DV_ZERO_M1_M1,   
    DV_ZERO_P1_M1,   DV_ZERO_M1_P1,

    DV_M1_P1_P1,       DV_P1_M1_M1,       
    DV_M1_M1_P1,       DV_P1_P1_M1,       
    DV_M1_M1_M1,       DV_P1_P1_P1,       
    DV_M1_P1_M1,       DV_P1_M1_P1
    }; 

    const amrex::GpuArray<int, constants::N_MICRO_STATES> bounce_Ydir={
    DV_ZERO_ZERO_ZERO, 
    DV_P1_ZERO_ZERO,   DV_M1_ZERO_ZERO,   
    DV_ZERO_M1_ZERO,   DV_ZERO_P1_ZERO,   
    DV_ZERO_ZERO_P1,   DV_ZERO_ZERO_M1,   
    
    DV_P1_M1_ZERO,   DV_M1_P1_ZERO,    
    DV_P1_P1_ZERO,   DV_M1_M1_ZERO,   
    DV_P1_ZERO_P1,   DV_M1_ZERO_M1,   
    DV_P1_ZERO_M1,   DV_M1_ZERO_P1,   
    DV_ZERO_M1_P1,   DV_ZERO_P1_M1,   
    DV_ZERO_M1_M1,   DV_ZERO_P1_P1, 
    
    DV_P1_M1_P1,       DV_M1_P1_M1,       
    DV_P1_P1_P1,       DV_M1_M1_M1,       
    DV_P1_P1_M1,       DV_M1_M1_P1,       
    DV_P1_M1_M1,       DV_M1_P1_P1
    };

const amrex::GpuArray<int, constants::N_MICRO_STATES> bounce_Zdir={
    DV_ZERO_ZERO_ZERO, 
    DV_P1_ZERO_ZERO,   DV_M1_ZERO_ZERO,   
    DV_ZERO_P1_ZERO,   DV_ZERO_M1_ZERO,   
    DV_ZERO_ZERO_M1,   DV_ZERO_ZERO_P1,

    DV_P1_P1_ZERO,   DV_M1_M1_ZERO,    
    DV_P1_M1_ZERO,   DV_M1_P1_ZERO,   
    DV_P1_ZERO_M1,   DV_M1_ZERO_P1,   
    DV_P1_ZERO_P1,   DV_M1_ZERO_M1,   
    DV_ZERO_P1_M1,   DV_ZERO_M1_P1,   
    DV_ZERO_P1_P1,   DV_ZERO_M1_M1, 
    
    DV_P1_P1_M1,       DV_M1_M1_P1,       
    DV_P1_M1_M1,       DV_M1_P1_P1,       
    DV_P1_M1_P1,       DV_M1_P1_M1,       
    DV_P1_P1_P1,       DV_M1_M1_M1
};


#endif
};

void check_stencil();
} // namespace lbm::stencil
#endif
