const std = @import("std");
const testing = std.testing;
const MultiPermute = @import("./MultiPermute.zig");

// Counting test

test "counts" {
    var start: usize = 0; // "runtime" for genuine slices
    _ = &start; // suppress warning/error

    const symbols3: [3]u8 = .{ 0, 0, 1 };
    var count = MultiPermute.countMultiset(symbols3[start..]);
    try testing.expectEqual(3, count);

    var symbols4: [4]u8 = .{ 0, 1, 1, 1 };
    count = MultiPermute.countMultiset(symbols4[start..]);
    try testing.expectEqual(4, count);

    symbols4 = .{ 1, 1, 2, 4 };
    count = MultiPermute.countMultiset(symbols4[start..]);
    try testing.expectEqual(12, count);

    const symbols6: [6]u8 = .{ 0, 2, 1, 3, 1, 1 };
    count = MultiPermute.countMultiset(symbols6[start..]);
    try testing.expectEqual(120, count);

    const symbols9: [9]u8 = .{ 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    count = MultiPermute.countMultiset(symbols9[start..]);
    try testing.expectEqual(362880, count);
}

// Sequence test

var visitor_i: usize = 0;
var visitor_item1 = [_]u8{0} ** MultiPermute.MaxN;
var visitor_item11 = [_]u8{0} ** MultiPermute.MaxN;

fn testseq_visitor(seq: []const u8) void {
    if (visitor_i == 1) {
        for (0..seq.len) |i| {
            visitor_item1[i] = seq[i];
        }
    } else if (visitor_i == 11) {
        for (0..seq.len) |i| {
            visitor_item11[i] = seq[i];
        }
    }
    visitor_i += 1;
}

test "sequence" {
    const symbols: [4]u8 = .{ 1, 1, 2, 4 };
    var start: usize = 0; // "runtime" for genuine slices
    _ = &start; // suppress warning/error

    MultiPermute.visitMultiset(symbols[start..], &testseq_visitor);

    const expected1: [4]u8 = .{ 2, 1, 1, 4 };
    try testing.expectEqual(expected1[0], visitor_item1[0]);
    try testing.expectEqual(expected1[1], visitor_item1[1]);
    try testing.expectEqual(expected1[2], visitor_item1[2]);
    try testing.expectEqual(expected1[3], visitor_item1[3]);

    const expected11: [4]u8 = .{ 1, 4, 2, 1 };
    try testing.expectEqual(expected11[0], visitor_item11[0]);
    try testing.expectEqual(expected11[1], visitor_item11[1]);
    try testing.expectEqual(expected11[2], visitor_item11[2]);
    try testing.expectEqual(expected11[3], visitor_item11[3]);
}

// TODO test full API with conserved prefix
// test "sequence" {
//     const N = 4; // the rightmost 4 symbols
//     const symbols = @Vector(16, u8){ 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 1, 1, 2, 4 };

//     MultiPermute.visitMultiset(symbols, N, &testseq_visitor);

//     // expecting the leftmost 12 to be conserved, and the remainder to follow a deterministic sequence
//     const expected1 = @Vector(16, u8){ 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 2, 1, 1, 4 };
//     try testing.expectEqual(expected1, visitor_item1);

//     const expected11 = @Vector(16, u8){ 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 1, 4, 2, 1 };
//     try testing.expectEqual(expected11, visitor_item11);
// }
