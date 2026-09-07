// Compile-time layout guards for the AoS peak types published as zero-copy
// structured numpy arrays by peaks_struct().
//
// The dtype overlays Peak1D / ChromatogramPeak / MobilityPeak1D directly, so it
// has to name the true byte offsets of position_ and intensity_. Both members
// are protected, which puts them out of reach of a plain offsetof from a binding
// translation unit. A derived probe reaches them: a derived class may access a
// protected member through its own type, and a member function body is a
// complete-class context, which is where offsetof needs the type to be complete.
//
// The probe declares no members of its own, so the peak stays the only
// data-holding class in the hierarchy. That keeps the probe standard-layout,
// which puts the peak subobject at offset 0 and makes the probe's offsets the
// peak's offsets; the assertion below checks that rather than assuming it.
//
// Marking one of the peak classes final would stop this header compiling. That
// is a build error with a confusing message rather than a silent problem, but
// it is the one edit that breaks the technique.

#pragma once

#include <cstddef>
#include <type_traits>

namespace pyopenms
{
    /// Byte offsets and ABI guards for a 1D peak type with protected position_/intensity_.
    ///
    /// Instantiating this template is what runs the checks; use its constants to
    /// build the dtype so the published offsets are the ones that were asserted.
    template <typename Peak>
    struct PeakLayout
    {
        struct Probe : Peak
        {
            static constexpr std::size_t position() { return offsetof(Probe, position_); }
            static constexpr std::size_t intensity() { return offsetof(Probe, intensity_); }
        };

        static constexpr std::size_t position_offset = Probe::position();
        static constexpr std::size_t intensity_offset = Probe::intensity();
        static constexpr std::size_t itemsize = sizeof(Peak);

        static_assert(std::is_standard_layout_v<Peak>,
            "peak must be standard-layout for zero-copy struct views (guarantees member order matches dtype)");
        // Standard-layout is what puts the peak subobject at offset 0, so the probe's
        // offsets are the peak's. sizeof(Probe) is deliberately not compared against
        // sizeof(Peak): an implementation may give an empty derived class tail padding
        // without moving any member, and a genuinely wrong offset is caught below.
        static_assert(std::is_standard_layout_v<Probe>,
            "layout probe must stay standard-layout, or its offsets do not describe the peak");
        static_assert(itemsize == 16,
            "peak must be 16 bytes for zero-copy structured array access");
        static_assert(std::is_same_v<typename Peak::CoordinateType, double>,
            "peak CoordinateType must be double (dtype assumes float64 for position)");
        static_assert(std::is_same_v<typename Peak::IntensityType, float>,
            "peak IntensityType must be float (dtype assumes float32 for intensity)");
        static_assert(position_offset == 0,
            "position_ must be the first member (dtype places the position field at offset 0)");
        static_assert(intensity_offset == 8,
            "intensity_ must directly follow the position with no padding (dtype places it at offset 8)");
    };
}
