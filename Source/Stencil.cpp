#include "Stencil.H"

namespace lbm::stencil {

void CheckStencil()
{
    const Stencil stencil;
    const auto& evs = stencil.evs;
    const auto& weight = stencil.weights;
    const auto& bounce_dirs = stencil.bounce_dirs;
    for (int q = 0; q < constants::n_micro_states; q++) {
        const int bounce_q = bounce_dirs[q];

        const auto& ev = evs[q];
        const auto& evb = evs[bounce_q];
        const int sum = abs(ev[0]) + abs(ev[1]) + abs(ev[2]);

        if ((ev[0] + evb[0] != 0) || (ev[1] + evb[1] != 0) ||
            (ev[2] + evb[2] != 0)) {
            amrex::Abort("Invalid bounce direction");
        }

        if (sum == 3 && abs(weight[q] - 1.0 / 216.0) > 1.0e-6) {
            amrex::Abort("Correct weight not set to 1/216");
        } else if (sum == 2 && abs(weight[q] - 1.0 / 54.0) > 1.0e-6) {
            amrex::Abort("Correct weight not set to 1/54");
        } else if (sum == 1 && abs(weight[q] - 2.0 / 27.0) > 1.0e-6) {
            amrex::Abort("Correct weight not set to 2/27");
        } else if (sum == 0 && abs(weight[q] - 8.0 / 27.0) > 1.0e-6) {
            amrex::Abort("Correct weight not set to 8/27");
        }
    }
}
} // namespace lbm::stencil
