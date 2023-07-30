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

test "sincos functions" {
    const Q_PI = 0x8000;
    const PI = std.math.pi;

    const FuncType: type = fn (u16) callconv(.C) SinCosExtern;
    const test_cases: [2]struct { FuncType, comptime_float, comptime_float } = .{
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
    const angle: f32 = castAngleIntToFloat(f32, @as(u16, x_folded));

    // Approximate positive cosine via Bhaskara's method
    //     => full range = [0, 1/2 pi]
    const cos = atoms.cosBhaskara(f32, angle);
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
    const x4th_bits = @as(composites.QuarterAngleAndHexadectrant, @bitCast(x_bits.angle_4th));

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

    fn cosBhaskara(comptime fN: type, x: fN) fN {
        // https://en.wikipedia.org/wiki/Bh%C4%81skara_I%27s_sine_approximation_formula#Equivalent_forms_of_the_formula
        const PI2 = comptime std.math.pi * std.math.pi;
        const xx = x * x;
        return (PI2 - mulPow2(2, xx)) / (PI2 + xx);
    }
};

fn SinCos(comptime fN: type) type {
    return struct {
        sin: fN,
        cos: fN,
    };
}

test "castAngleIntToFloat" {
    const Q_EIGHTH_PI = 0x1000;
    const EIGHTH_PI = std.math.pi / 8.0;
    const test_cases: [16]struct { f32, u16 } = .{
        .{ 0, 0 },
        .{ EIGHTH_PI, Q_EIGHTH_PI },
        .{ 2 * EIGHTH_PI, 2 * Q_EIGHTH_PI },
        .{ 3 * EIGHTH_PI, 3 * Q_EIGHTH_PI },
        .{ 4 * EIGHTH_PI, 4 * Q_EIGHTH_PI },
        .{ 5 * EIGHTH_PI, 5 * Q_EIGHTH_PI },
        .{ 6 * EIGHTH_PI, 6 * Q_EIGHTH_PI },
        .{ 7 * EIGHTH_PI, 7 * Q_EIGHTH_PI },
        .{ 8 * EIGHTH_PI, 8 * Q_EIGHTH_PI },
        .{ 9 * EIGHTH_PI, 9 * Q_EIGHTH_PI },
        .{ 10 * EIGHTH_PI, 10 * Q_EIGHTH_PI },
        .{ 11 * EIGHTH_PI, 11 * Q_EIGHTH_PI },
        .{ 12 * EIGHTH_PI, 12 * Q_EIGHTH_PI },
        .{ 13 * EIGHTH_PI, 13 * Q_EIGHTH_PI },
        .{ 14 * EIGHTH_PI, 14 * Q_EIGHTH_PI },
        .{ 15 * EIGHTH_PI, 15 * Q_EIGHTH_PI },
    };

    for (test_cases) |test_case| {
        const result = castAngleIntToFloat(f32, test_case[1]);
        const expt_result = test_case[0];

        try std.testing.expectApproxEqRel(expt_result, result, 1e-6);
    }
}

/// Assumes the full range of @TypeOf(angle) linearly spans [0, 2*pi].
fn castAngleIntToFloat(comptime f: type, angle: anytype) f {
    const IntType = @TypeOf(angle);
    const scale_to_radians = _2PI / @as(
        comptime_float,
        @floatFromInt(1 << @bitSizeOf(IntType)),
    );

    var angle_f: f = @floatFromInt(angle);
    return angle_f * scale_to_radians;
}

test "mulPow2" {
    try std.testing.expectEqual(@as(f32, 2.0), mulPow2(1, @as(f32, 1.0)));
    try std.testing.expectEqual(@as(f32, 0.5), mulPow2(-1, @as(f32, 1.0)));

    const shift_to_pos_zero = mulPow2(-5, std.math.floatMin(f32));
    try std.testing.expect(mulPow2(0, std.math.floatMin(f32)) > 0);
    try std.testing.expectEqual(@as(f32, 0), shift_to_pos_zero);
    try std.testing.expectEqual(false, std.math.signbit(shift_to_pos_zero));

    const shift_to_neg_zero = mulPow2(-5, -std.math.floatMin(f32));
    try std.testing.expect(mulPow2(0, -std.math.floatMin(f32)) < 0);
    try std.testing.expectEqual(@as(f32, -0), shift_to_neg_zero);
    try std.testing.expectEqual(true, std.math.signbit(shift_to_neg_zero));
}

fn mulPow2(
    comptime pow2: comptime_int,
    x: anytype,
) @TypeOf(x) {
    const fN = @TypeOf(x);
    comptime std.debug.assert(@typeInfo(fN) == .Float);

    // Short-circuit no-ops
    if (comptime pow2 == 0) return x;

    var x_bits: FloatDecompSign(fN) = @bitCast(x);
    const exp_1 = comptime @as(FloatDecompSign(fN), @bitCast(
        FloatDecomp(fN){ .mantissa = 0, .exponent = 1, .sign_neg = false },
    ));

    if (pow2 > 0) {
        // TODO handle overflow case
        x_bits.exp_man += pow2 * exp_1.exp_man;
    } else {
        // Saturating subtraction: if power-of-two denominator is too large,
        // the result will set all matissa bits to zero -> float will be zero.
        x_bits.exp_man -|= (-pow2) * exp_1.exp_man;
    }
    return @bitCast(x_bits);
}

test "floating point decomp" {
    inline for (.{ f16, f32, f64 }) |f| {
        comptime try std.testing.expectEqual(
            @bitSizeOf(f),
            @bitSizeOf(FloatDecomp(f)),
        );

        const test_cases: [3]struct { f, FloatDecomp(f) } = .{
            .{ 0.0, .{
                .sign_neg = false,
                .exponent = 0,
                .mantissa = 0,
            } },
            .{ 1.0, .{
                .sign_neg = false,
                .exponent = (1 << (std.math.floatExponentBits(f) - 1)) - 1,
                .mantissa = 0,
            } },
            .{ -0.375, .{
                .sign_neg = true,
                .exponent = (1 << (std.math.floatExponentBits(f) - 1)) - 3,
                .mantissa = 1 << (std.math.floatMantissaBits(f) - 1),
            } },
        };
        inline for (test_cases) |test_case| {
            const result: FloatDecomp(f) = @bitCast(test_case[0]);
            const expt_result = test_case[1];

            try std.testing.expectEqual(expt_result, result);
        }
    }
}

fn FloatDecomp(comptime f: type) type {
    return packed struct {
        mantissa: std.meta.Int(std.builtin.Signedness.unsigned, std.math.floatMantissaBits(f)),
        exponent: std.meta.Int(std.builtin.Signedness.unsigned, std.math.floatExponentBits(f)),
        sign_neg: bool,
    };
}

fn FloatDecompSign(comptime f: type) type {
    return packed struct {
        exp_man: std.meta.Int(
            std.builtin.Signedness.unsigned,
            std.math.floatExponentBits(f) + std.math.floatMantissaBits(f),
        ),
        sign_neg: bool,
    };
}

const Q_8TH_PI = 0x1000;
const Q_4TH_PI = 0x2000;
const Q_HALF_PI = 0x4000;

const SQRT_HALF = std.math.sqrt1_2;
const _2PI = 2 * std.math.pi;
