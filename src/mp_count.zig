const builtin = @import("builtin");
const std = @import("std");
const MultiPermute = @import("./MultiPermute.zig");

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

    const count = MultiPermute.countMultiset(symbols, i);
    std.debug.print("{d}\n", .{count});
}
