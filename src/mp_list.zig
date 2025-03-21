const builtin = @import("builtin");
const std = @import("std");
const MultiPermute = @import("./MultiPermute.zig");

var n: usize = 0;
const stdout = std.io.getStdOut().writer();

pub fn main() !void {
    var symbols = [_]u8{0} ** MultiPermute.MaxN;

    var args = std.process.args(); // not on windows
    var i: usize = 0;
    _ = args.next();

    readargs: while (args.next()) |arg| {
        for (arg) |c| {
            if (i < MultiPermute.MaxN) symbols[i] = c;
            i += 1;
            if (i > MultiPermute.MaxN) break :readargs;
        }
    }

    if (i > MultiPermute.MaxN) {
        std.debug.print("Max {d} symbols.\n", .{MultiPermute.MaxN});
        return;
    }

    n = i;

    MultiPermute.visitMultiset(symbols[0..i], &logMultiset);
}

var symbolIO = [_]u8{0} ** 16;

fn logMultiset(seq: []const u8) void {
    stdout.print("{s}\n", .{seq[0..]}) catch unreachable;
}
