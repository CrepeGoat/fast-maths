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
    const sample_count = 64;
    const coord_low = -10.0;
    const coord_high = 10.0;

    const tolerance = 9607;

    var max: u16 = 0;
    for (0..sample_count) |i| for (0..sample_count) |j| {
        const n: comptime_float = @floatFromInt(sample_count);
        const i_float: f32 = @floatFromInt(i);
        const j_float: f32 = @floatFromInt(j);
        const x: f32 = ((coord_high - coord_low) / n) * i_float - coord_low;
        const y: f32 = ((coord_high - coord_low) / n) * j_float - coord_low;

        const result = atan2Cordic3Poly1(y, x);
        const expt_result: u16 = @intFromFloat(std.math.atan2(f32, y, x) * (1 << @bitSizeOf(u16)) / (2 * std.math.pi));

        const epsilon: i16 = @bitCast(expt_result -% result);
        max = @max(max, std.math.absCast(epsilon));
        try std.testing.expect(std.math.absCast(epsilon) <= tolerance);
    };
    std.debug.print("max: {d}\n", .{max});
}

pub export fn atan2Rational2(y: f32, x: f32) u16 {
    var x_mut = x;
    var y_mut = y;

    const x_abs = @fabs(x);
    const y_abs = @fabs(y);

    const is_refl_pi_4th = (x_abs > y_abs);
    if (is_refl_pi_4th) {
        // We save 1 MOV instruction by using a -temp here,
        // which performs a NEG and stores the result in x's register.
        const tmp = y_mut;
        y_mut = x_mut;
        x_mut = -tmp;
    }

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

    const y2 = y_mut * y_mut;
    const x2 = x_mut * x_mut;

    const num = x_mut * (A * y2 + B * x2);
    const den = y_mut * (y2 + C * x2);
    if (den == 0.0) {
        // If the math is correct here, we don't have to check for NaNs or infinities.
        // Just the single case of denominator == 0.0. Otherwise we get back a valid number.
        // Also, there should be 4 more instructions from this if statement, and 2-4 more cycles.
        return 0;
    }

    // Do a FLOOR and MFC1 instruction combo here! That should add 6 cycles and 2 instructions.
    var angle: u16 = @bitCast(@as(i16, @intFromFloat(@floor(num / den))));

    // This sequence of conditions should be 6 instructions and 4-6 cycles total.
    if (is_refl_pi_4th) {
        angle +%= Q_HALF_PI;
    }
    if (y_mut < 0.0) {
        // Using some clever facts of how we swapped around y earlier!
        // This lets us do two rotational checks in one fell swoop.
        angle += Q_PI;
    }
    // Conversion to s16 adds 2 instructions and 2 cycles.
    // Return instruction is 1 more instruction and 1 more cycle.
    return angle;
}

pub export fn atan2Cordic3Poly1(y: f32, x: f32) u16 {
    var x_mut = x;
    var y_mut = y;

    const is_rot_180 = y_mut < 0;
    if (is_rot_180) {
        // Rotate -180
        const tmp = x_mut;
        x_mut = -y_mut;
        y_mut = -tmp;
    }
    const is_rot_90 = x_mut < 0;
    if (is_rot_90) {
        // Rotate -90
        const tmp = x_mut;
        x_mut = y_mut;
        y_mut = -tmp;
    }
    const is_rot_45 = x_mut < y_mut;
    if (is_rot_45) {
        // Rotate -45
        const x_tmp = x_mut;
        const y_tmp = y_mut;
        x_mut = SQRT_HALF * (x_tmp + y_tmp);
        y_mut = SQRT_HALF * (y_tmp - x_tmp);
    }

    if (x_mut == 0) {
        return 0;
    }

    const radians_to_uint = @as(comptime_float, 1 << @bitSizeOf(u16)) / _2PI;
    const slope_adj = 0.5 * (1 + SQRT_HALF);
    var result: u16 = @intFromFloat((comptime radians_to_uint * slope_adj) * y_mut / x_mut);

    if (is_rot_45) {
        result += Q_4TH_PI;
    }
    if (is_rot_90) {
        result += Q_HALF_PI;
    }
    if (is_rot_180) {
        result += Q_PI;
    }

    return result;
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
