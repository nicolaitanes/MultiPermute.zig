const builtin = @import("builtin");
const std = @import("std");
const MultiPermute = @import("./MultiPermute.zig");

var n: usize = 0;
const stdout = std.io.getStdOut().writer();

pub fn main() !void {
    var symbols: @Vector(16, u8) = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    var args = std.process.args(); // not on windows
    var i: usize = 0;
    _ = args.next();

    readargs: while (args.next()) |arg| {
        for (arg) |c| {
            if (i < 16) symbols[15 - i] = c;
            i += 1;
            if (i > 16) break :readargs;
        }
    }

    if (i > 16) {
        std.debug.print("Max 16 symbols.\n", .{});
        return;
    }

    n = i;

    MultiPermute.visitMultiset(symbols, i, &logMultiset);
}

var symbolIO = [_]u8{0} ** 16;

fn logMultiset(seq: @Vector(16, u8)) void {
    symbolIO[0..16].* = seq;
    stdout.print("{s}\n", .{symbolIO[16 - n ..]}) catch unreachable;
}
