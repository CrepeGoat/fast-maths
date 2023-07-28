// zig version 0.11.0-dev.4059+17255bed4
const std = @import("std");

test "fibonacci numbers 0 to 10" {
    const expt_values = [_]u64{ 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55 };

    for (expt_values, 0..expt_values.len) |expt_value, i| {
        const result = fibonacci(u64, make_int_context(u64), i);
        try std.testing.expectEqual(expt_value, result);
    }
}

test "fibonacci number one million" {
    const context = struct {
        const mod = 1000007;
        pub fn make_zero() u64 {
            return 0;
        }
        pub fn make_one() u64 {
            return 1;
        }
        pub fn add(x1: u64, x2: u64) u64 {
            return (x1 + x2) % mod;
        }
        pub fn mul(x1: u64, x2: u64) u64 {
            return (x1 * x2) % mod;
        }
    };

    const N = 1000000;
    try std.testing.expectEqual(
        fibonacci(u64, context, @as(u64, N)),
        context.add(
            fibonacci(u64, context, @as(u64, N - 1)),
            fibonacci(u64, context, @as(u64, N - 2)),
        ),
    );
}

pub fn fibonacci(
    comptime NumType: type,
    comptime context: anytype,
    n: anytype,
) NumType {
    const NType = @TypeOf(n);

    var buffer: [3]NumType = .{
        context.make_one(),
        context.make_zero(),
        context.make_one(),
    };
    var buffer_extra: @TypeOf(buffer) = undefined;

    var buf1 = &buffer;
    var buf2 = &buffer_extra;

    var index_bit = @bitSizeOf(NType) - @clz(n);
    while (index_bit > 0) {
        index_bit -= 1;
        const bit = (n >> @intCast(index_bit)) & 1 != 0;
        const index_offset: u8 = if (bit) 1 else 0;

        for (0..2) |i| {
            const j = i + index_offset;

            buf2[i] = context.add(
                context.mul(buf1[j / 2], buf1[(j + 1) / 2]),
                context.mul(buf1[j / 2 + 1], buf1[(j + 1) / 2 + 1]),
            );
        }
        buf2[2] = context.add(buf2[0], buf2[1]);

        std.mem.swap(*[3]NumType, &buf1, &buf2);
    }

    return buf1[1];
}

pub fn make_int_context(comptime Int: type) type {
    comptime std.debug.assert(@typeInfo(Int) == .Int);

    return struct {
        pub fn make_zero() Int {
            return 0;
        }
        pub fn make_one() Int {
            return 1;
        }
        pub fn add(x1: Int, x2: Int) Int {
            return x1 + x2;
        }
        pub fn mul(x1: Int, x2: Int) Int {
            return x1 * x2;
        }
    };
}
