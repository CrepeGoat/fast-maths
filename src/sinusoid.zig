/// zig version 0.11.0-dev.4059+17255bed4
///
/// This is a collection of approximate sine & cosine functions of minimal cycle
/// duration.
///
/// This was originally inspired by Kaze Emanuar's work on optimizing SM64 (see
/// https://youtu.be/xFKFoGiGlXQ).
/// As Kaze's work is about execution on the N64, the comments on optimization
/// in this file are based on the assembly output when compiling for the
/// `mips-linux-gnu` target.
const std = @import("std");

comptime {
    @setFloatMode(std.builtin.FloatMode.Optimized);
}

test "arctan2 functions" {
    const u16_count = 1 << @bitSizeOf(u16);
    const sample_count = 256;
    const radii = [_]f32{ 1e-7, 1, 1e7 };

    const FuncType: type = fn (f32, f32) callconv(.C) u16;
    const test_cases = [_]struct { FuncType, comptime_int }{
        .{ atan2Cordic2Rational3, 2 },
        .{ atan2Cordic2Poly1, 500 },
    };

    inline for (test_cases) |test_case| {
        const func = test_case[0];
        const tolerance = test_case[1];

        var max_eps: u16 = 0;
        for (0..sample_count) |i| {
            const expt_result: u16 = @as(u16, @intCast(i)) * (u16_count / sample_count);
            const expt_angle = @as(f32, @floatFromInt(expt_result)) * _2PI / @as(f32, u16_count);

            for (radii) |radius| {
                const x: f32 = std.math.cos(expt_angle) * radius;
                const y: f32 = std.math.sin(expt_angle) * radius;
                const result: i16 = @bitCast(func(y, x));

                const epsilon = @as(i16, @bitCast(expt_result)) -% result;
                max_eps = @max(max_eps, std.math.absCast(epsilon));
                try std.testing.expect(std.math.absCast(epsilon) <= tolerance);
            }
        }

        // Corner-case: <x, y> = <0, 0>
        try std.testing.expectEqual(@as(u16, 0), func(0.0, 0.0));
    }
}

pub export fn atan2Cordic2Rational3(y: f32, x: f32) u16 {
    const atan2At0 = struct {
        fn call(_y: f32, _x: f32) u16 {
            return atoms.atan2Rational3At0Mod180(f32, _y, _x);
        }
    }.call;

    return composites.atan2Cordic2(f32, atan2At0, y, x) catch 0.0;
}

pub export fn atan2Cordic2Poly1(y: f32, x: f32) u16 {
    const atan2At0 = struct {
        fn call(_y: f32, _x: f32) u16 {
            return atoms.atan2Poly1At0Mod180(f32, _y, _x);
        }
    }.call;

    return composites.atan2Cordic2(f32, atan2At0, y, x) catch 0.0;
}

test "sincos functions" {
    const PI = std.math.pi;

    const FuncType: type = fn (u16) callconv(.C) SinCosExtern;
    const test_cases = [_]struct { FuncType, comptime_float, comptime_float }{
        .{ sincosBhaskara, 3e-3, 2e-7 },
        .{ sincosPoly3Cordic1, 3e-3, 6e-5 },
    };
    inline for (test_cases) |test_case| {
        const sincos = test_case[0];
        const sincosEpsilon = test_case[1];
        const pythagoreanEpsilon = test_case[2];

        const sample_count: u16 = 256;
        for (0..sample_count) |i| {
            const angle_uint: u16 = @as(u16, @truncate(i)) * 2 * (Q_PI / sample_count);
            const angle = @as(f32, @floatFromInt(i)) * (2 * PI / @as(
                comptime_float,
                @floatFromInt(sample_count),
            ));

            const result = sincos(angle_uint);
            const expt_sin = std.math.sin(angle);
            const expt_cos = std.math.cos(angle);

            // Confirm individual sine/cosine values are nearly accurate
            try std.testing.expectApproxEqAbs(expt_sin, result.sin, sincosEpsilon);
            try std.testing.expectApproxEqAbs(expt_cos, result.cos, sincosEpsilon);

            // Confirm sin^2 + cos^2 is near 1
            try std.testing.expectApproxEqAbs(
                @as(f32, 1),
                std.math.pow(f32, result.sin, 2) + std.math.pow(f32, result.cos, 2),
                pythagoreanEpsilon,
            );
        }
    }
}

pub export fn sincosBhaskara(x: u16) SinCosExtern {
    const x_bits = @as(composites.AngleAndQuadrant, @bitCast(x));
    const cos_lt0 = x_bits.is_even_half != x_bits.is_even_4th;
    const sin_lt0 = x_bits.is_even_half;

    // OPTIMIZATION: another alternative is to condition on cos_lt0:
    //
    // var x_folded: i16 = @bitCast(x);
    // if (cos_lt0) {
    //     x_folded = std.math.minInt(i16) -% x_folded;
    // }
    //
    // Using the below expression instead saves 1 MIPS instruction.
    // ---
    // Strip the most significant 2 bits of the fixed-point angle
    //     [0, 4/2 pi] -> [0, 1/2 pi]
    //     => full range = [0, 1/2 pi]
    var x_folded: u15 = x_bits.angle_4th;
    // Conditionally reflect angle about theta = 1/4 pi
    //     [0, 1/2 pi] -> [1/2 pi, 0]
    //     => full range = [0, 1/2 pi]
    if (x_bits.is_even_4th) {
        x_folded = Q_HALF_PI -% x_folded;
    }

    // Approximate positive cosine via Bhaskara's method
    //     => full range = [0, 1/2 pi]
    const cos = atoms.cosBhaskara(f32, @as(u16, x_folded));
    // Calculate positive sine from cosine
    const sin = composites.posSinFromCos(f32, cos);

    // Conditionally adjust signs of the result based on the angle's quadrant
    //     [0, 1/2 pi] -> [0, 4/2 pi]
    //     => full range = [0, 2 pi]
    return .{
        .sin = if (sin_lt0) -sin else sin,
        .cos = if (cos_lt0) -cos else cos,
    };
}

pub export fn sincosPoly3Cordic1(x: u16) SinCosExtern {
    const x_bits = @as(composites.AngleAndQuadrant, @bitCast(x));
    const x4th_bits = @as(
        composites.QuarterAngleAndHexadectrant,
        @bitCast(x_bits.angle_4th),
    );

    const cos_lt0 = x_bits.is_even_half != x_bits.is_even_4th;
    const sin_lt0 = x_bits.is_even_half;
    const is_refl_pi_half = x_bits.is_even_4th;
    const is_refl_pi_4th = x4th_bits.is_even_16th;
    const is_rot_pi_4th = x4th_bits.is_even_8th != x4th_bits.is_even_16th;
    const is_refl = is_refl_pi_half != is_refl_pi_4th;

    // Strip the most significant 4 bits of the fixed-point angle
    //     [0, 16/8 pi] -> [0, 1/8 pi]
    //     => full range = [0, 1/8 pi]
    var x_folded: u13 = x4th_bits.angle_16th;
    // Conditionally reflect angle about theta = 1/16 pi
    //     [0, 1/8 pi] -> [1/8 pi, 0]
    //     => full range = [0, 1/8 pi]
    if (is_refl_pi_4th) {
        x_folded = Q_8TH_PI - x_folded;
    }

    // Approximate sin & cos via polynomials
    var result = atoms.sincosPoly3ApproxComplement(f32, @as(u16, x_folded));

    // Conditionally rotate the result 1/4 pi radians
    //     [0, 1/8 pi] -> [2/8 pi, 3/8 pi]
    //     => full range = [0, 1/8 pi] && [2/8 pi, 3/8 pi]
    if (is_rot_pi_4th) {
        const sin_tmp = SQRT_HALF * (result.cos + result.sin);
        const cos_tmp = SQRT_HALF * (result.cos - result.sin); // sin < cos -> always pos
        result.sin = sin_tmp;
        result.cos = cos_tmp;
    }
    // Conditionally reflect the result about theta = 1/4 pi radians
    //     [0, 1/8 pi] -> [4/8 pi, 3/8 pi]
    //     [2/8 pi, 3/8 pi] -> [2/8 pi, 1/8 pi]
    //     => full range = [0, 4/8 pi]
    if (is_refl) {
        std.mem.swap(f32, &result.sin, &result.cos);
    }
    // Conditionally adjust signs of the result based on the angle's quadrant
    //     [0, 4/8 pi] -> [0, 16/8 pi]
    //     => full range = [0, 2 pi]
    return .{
        .sin = if (sin_lt0) -result.sin else result.sin,
        .cos = if (cos_lt0) -result.cos else result.cos,
    };
}

pub const SinCosExtern = extern struct {
    sin: f32,
    cos: f32,
};

const composites = struct {
    fn atan2Cordic2(
        comptime Float: type,
        comptime atan2At0: fn (Float, Float) u16,
        y: Float,
        x: Float,
    ) error{AngleAtOrigin}!u16 {
        var x_mut = x;
        var y_mut = y;

        // Rotate input
        const is_refl_45 = @fabs(x_mut) < @fabs(y_mut);
        if (is_refl_45) {
            std.mem.swap(Float, &x_mut, &y_mut);
        }
        std.debug.assert(@fabs(x_mut) >= @fabs(y_mut));

        if (x_mut == 0) {
            return error.AngleAtOrigin;
        }

        const is_rot_180 = x_mut < 0;
        // tan(x) = tan(x + 180) -> no need to rotate 180

        // Approximate atan about tan = 0
        var result: u16 = atan2At0(y_mut, x_mut);

        // Un-rotate output
        if (is_rot_180) {
            result +%= Q_PI;
        }
        if (is_refl_45) {
            result = Q_HALF_PI -% result;
        }

        return result;
    }

    const AngleAndQuadrant = packed struct {
        angle_4th: u14,
        is_even_4th: bool,
        is_even_half: bool,
    };

    const QuarterAngleAndHexadectrant = packed struct {
        angle_16th: u12,
        is_even_16th: bool,
        is_even_8th: bool,
    };

    fn posSinFromCos(comptime fN: type, cos: fN) fN {
        return std.math.sqrt(1 - cos * cos);
    }
};

const atoms = struct {
    fn atan2Poly1At0Mod180(comptime fN: type, y: fN, x: fN) u16 {
        const radians_to_uint = @as(comptime_float, 1 << @bitSizeOf(u16)) / _2PI;
        const slope_adj = 5.0 / 6.0;
        return @bitCast(@as(
            i16,
            @intFromFloat((comptime radians_to_uint * slope_adj) * (y / x)),
        ));
    }

    fn atan2Rational3At0Mod180(comptime fN: type, y: fN, x: fN) u16 {
        // Math section follows.
        // In total, we perform 7 multiplications and 2 additions.
        // That should be 41 cycles total.
        // Including the division (29 cycles), that's 70 cycles of math ops.

        // These constants chosen for a minimax approx subject to constraints
        // that f(1) = 8192, i.e. the s16 angle corresponding to pi/4.
        // Results are valid for magnitudes of (x, y) between
        // 1.0 x 10^-13 and 1.0 x 10^11.
        const A = 10420.2706879;
        const B = 1956.92683519;
        const C = 0.510888369518;

        const y2 = y * y;
        const x2 = x * x;

        const num = y * (A * x2 + B * y2);
        const den = x * (x2 + C * y2);
        std.debug.assert(den != 0.0);

        // Do a FLOOR and MFC1 instruction combo here! That should add 6 cycles and 2 instructions.
        return @bitCast(@as(i16, @intFromFloat(@floor(num / den))));
    }

    fn sincosPoly3ApproxComplement(comptime fN: type, x: anytype) SinCos(fN) {
        const IntType = @TypeOf(x);
        comptime std.debug.assert(@typeInfo(IntType) == .Int);
        // Assumes the full range of @TypeOf(x) linearly spans [0, 2*pi].
        const scale_to_radians = _2PI / @as(
            comptime_float,
            @floatFromInt(1 << @bitSizeOf(IntType)),
        );

        const x_float: fN = @floatFromInt(x);

        // OPTIMIZATION: the original polynomials can be written as
        //
        //  sin(x) = x (1 - 1/8 x^2)
        //  cos(x) = 1 - 1/2 x^2
        //
        // By storing the integer angle in half-radians (1/2 x) instead of
        // radians (x), and then storing its square (1/4 x^2), we can compute
        // the same polynomials as
        //
        //  sin(x) = (1/2 x) (2 - (1/2 x)^2)
        //  cos(x) = 1 - (1/2 x)^2 - (1/2 x)^2
        //
        // Overall this uses 2 fewer multiplication instructions, in exchange
        // for 1 extra subtraction in the cosine calculation.
        const x_half = (comptime 0.5 * scale_to_radians) * x_float;
        const xx_4th = x_half * x_half;

        return .{
            .sin = x_half * (2 - xx_4th),
            .cos = 1 - xx_4th - xx_4th,
        };
    }

    fn cosBhaskara(comptime fN: type, x: anytype) fN {
        const IntType = @TypeOf(x);
        comptime std.debug.assert(@typeInfo(IntType) == .Int);

        const x_float: fN = @floatFromInt(x);

        // https://en.wikipedia.org/wiki/Bh%C4%81skara_I%27s_sine_approximation_formula#Equivalent_forms_of_the_formula
        // ---
        // OPTIMIZATION: the original rational expression can be written as:
        //
        //  f(x) = (PI^2 - 4x^2) / (PI^2 + x^2)
        //
        // where x is in radians. This can be refactored into a function that
        // operates on angles in arbitrary units u := s * radians:
        //
        //  f(x) = (s^2 * PI^2 - 4(s*x)^2) / (s^2 * PI^2 + (s*x)^2)
        //
        // By choosing s = (2^16) / (2PI), we can skip the scaling
        // multiplication performed after the int-to-float conversion.
        const U_PI = comptime @as(
            comptime_float,
            @floatFromInt(1 << @bitSizeOf(IntType)),
        ) / 2.0;
        const U_PI2 = comptime U_PI * U_PI;

        const xx = x_float * x_float;
        const _2xx = xx + xx;
        const _4xx = _2xx + _2xx;
        return (U_PI2 - _4xx) / (U_PI2 + xx);
    }
};

fn SinCos(comptime fN: type) type {
    return struct {
        sin: fN,
        cos: fN,
    };
}

const Q_8TH_PI = 0x1000;
const Q_4TH_PI = 0x2000;
const Q_HALF_PI = 0x4000;
const Q_PI = 0x8000;

const SQRT_HALF = std.math.sqrt1_2;
const _2PI = 2 * std.math.pi;
