const builtin = @import("builtin");
const std = @import("std");
const MultiPermute = @import("./MultiPermute.zig");

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

    const count = MultiPermute.countMultiset(symbols[0..i]);
    std.debug.print("{d}\n", .{count});
}
