const builtin = @import("builtin");
const std = @import("std");
const MultiPermute = @import("./MultiPermute.zig");

var n: usize = 0;
const stdout = std.io.getStdOut().writer();

pub fn main() !void {
    var letters: @Vector(16, u8) = .{ 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32 };

    var args = std.process.args(); // not on windows
    var i: usize = 0;
    _ = args.next();

    readargs: while (args.next()) |arg| {
        for (arg) |c| {
            if (i < 16) letters[15 - i] = c;
            i += 1;
            if (i > 16) break :readargs;
        }
    }

    if (i > 16) {
        std.debug.print("Max 16 letters.\n", .{});
        return;
    }

    n = i;
    letterIO[16] = 0;

    MultiPermute.visitMultiset(letters, i, &logMultiset);
}

var letterIO = [_]u8{0} ** 17;

fn logMultiset(seq: @Vector(16, u8)) void {
    letterIO[0..16].* = seq;
    stdout.print("{s}\n", .{letterIO[16 - n ..]}) catch unreachable;
}
