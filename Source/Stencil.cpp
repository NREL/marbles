#include "Stencil.H"

namespace lbm::stencil {

void check_stencil()
{
    const Stencil stencil;
    const auto& evs = stencil.evs;
    const auto& weight = stencil.weights;
    const auto& bounce_dirs = stencil.bounce_dirs;
    for (int q = 0; q < constants::N_MICRO_STATES; q++) {
        const int bounce_q = bounce_dirs[q];

        const auto& ev = evs[q];
        const auto& evb = evs[bounce_q];
        const int sum = std::abs(ev[0]) + std::abs(ev[1]) + std::abs(ev[2]);

        if ((ev[0] + evb[0] != 0) || (ev[1] + evb[1] != 0) ||
            (ev[2] + evb[2] != 0)) {
            amrex::Abort("Invalid bounce direction");
        }

        if ((sum == 3) &&
            (std::abs(weight[q] - 1.0 / 216.0) > constants::SMALL_NUM)) {
            amrex::Abort("Correct weight not set to 1/216");
        } else if (
            (sum == 2) &&
            (std::abs(weight[q] - 1.0 / 54.0) > constants::SMALL_NUM)) {
            amrex::Abort("Correct weight not set to 1/54");
        } else if (
            (sum == 1) &&
            (std::abs(weight[q] - 2.0 / 27.0) > constants::SMALL_NUM)) {
            amrex::Abort("Correct weight not set to 2/27");
        } else if (
            (sum == 0) &&
            (std::abs(weight[q] - 8.0 / 27.0) > constants::SMALL_NUM)) {
            amrex::Abort("Correct weight not set to 8/27");
        }
    }
}
} // namespace lbm::stencil
