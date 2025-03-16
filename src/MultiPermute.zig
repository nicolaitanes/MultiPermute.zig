// This module encodes functions to generate the permutations of a multiset
// following this algorithm:

// Algorithm 1
// Visits the permutations of multiset E. The permutations are stored
// in a singly-linked list pointed to by head pointer h. Each node in the linked
// list has a value field v and a next field n. The init(E) call creates a
// singly-linked list storing the elements of E in non-increasing order with h, i,
// and j pointing to its first, second-last, and last nodes, respectively. The
// null pointer is given by φ. Note: If E is empty, then init(E) should exit.
// Also, if E contains only one element, then init(E) does not need to provide a
// value for i.

// A linked list must be used instead of an array because, as long as we keep a
// reference to the nodes at the removal and insertion points, arbitrary shifts
// (rotations of sub-lists) can be performed in constant time; this is what h, i,
// and j are for.

// [h, i, j] ← init(E)
// visit(h)
// while j.n ≠ φ or j.v < h.v do
//     if j.n ≠ φ and i.v ≥ j.n.v then
//         s←j
//     else
//         s←i
//     end if
//     t←s.n
//     s.n ← t.n
//     t.n ← h
//     if t.v < h.v then
//         i←t
//     end if
//     j←i.n
//     h←t
//     visit(h)
// end while

// Adapted from https://github.com/ekg/multipermute
// and optimized for fast WebAssembly, with a limit 3-16 elements in the multiset.

// The MIT License (MIT)

// Copyright (c) 2025 Christopher Nicolai

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

const std = @import("std");
const testing = std.testing;

pub const MaxN = 16;
const zero128 = @Vector(4, u32){ 0, 0, 0, 0 };

pub fn MultiPermute(comptime N: u32) type {
    return MultiPermuteWithVisitor(N, struct {});
}

// For applications where reading each permutation (traversing a linked list)
// is a significant fraction of runtime, provide a Visitor that can matchQuad(fourBytes, index in 0..4) --
// if not matched, list traversal is abandoned and we move on to the next permutation.
pub fn MultiPermuteWithVisitor(comptime N: u32, comptime Visitor: type) type {
    if (N > MaxN) @compileError("Max 16 symbols.");
    return struct {
        /// Visitor: { .matchQuad(x, i) bool, .match(x) bool }
        visitor: Visitor,
        nodes: [N]Node,
        permu: [N]u8,
        remap: [N]u8,
        N: u32 = N,
        h: ?*Node = null,
        i: ?*Node = null,
        j: ?*Node = null,
        input: @Vector(16, u8) = @bitCast(zero128),

        const Self = @This();
        const Node = struct { v: u32, r: u32, n: ?*Node };

        pub inline fn initVisitor(visitor: Visitor) Self {
            return Self{ .visitor = visitor, .nodes = [_]Node{Node{ .v = 0, .r = 0, .n = null }} ** N, .permu = [_]u8{0} ** N, .remap = [_]u8{0} ** N };
        }

        pub inline fn init() Self {
            return initVisitor(.{});
        }

        pub inline fn reset(self: *Self) void {
            self.h = null;
        }

        // Inputs symbols for permutation. The first (N-16) symbols will be untouched in the output;
        // the remaining N symbols will be permuted.
        // adapted from ms2permutation
        pub fn initLetters(self: *Self, letters: @Vector(16, u8)) void {
            self.input = letters;
            // ms2multiplicities:
            var counters = [_]u8{0} ** 256;
            for (0..N) |i| {
                const e = letters[MaxN - i - 1];
                counters[e] += 1;
            }
            var mults = [_]u8{0} ** N;
            var nMult: usize = 0;
            for (0..N) |i| {
                const e = letters[MaxN - i - 1];
                const n = counters[e];
                if (n == 0) continue;
                self.remap[nMult] = e;
                mults[nMult] = n;
                counters[e] = 0;
                nMult += 1;
            }
            // multiplicity2sorted:
            var j: usize = 0;
            for (0..N) |k| {
                var m = mults[k];
                while (m > 0) : (m -= 1) {
                    self.permu[j] = @intCast(k);
                    j += 1;
                }
            }
            self.reset();
        }

        pub fn initNodes(self: *Self) void {
            if (self.permu[N - 1] == 0)
                return self.initLetters(self.input);

            for (0..N) |i| {
                self.nodes[i].v = self.permu[N - i - 1];
                const v: u32 = @bitCast(self.nodes[i].v);
                self.nodes[i].r = self.remap[v];
                self.nodes[i].n = if (i < N - 1) &self.nodes[i + 1] else null;
            }
            self.h = &self.nodes[0];
            self.i = &self.nodes[N - 2];
            self.j = &self.nodes[N - 1];
        }

        pub inline fn nextFrom(self: *Self, h: ?*Node, i: ?*Node, j: ?*Node) [3]?*Node {
            var s = i.?;

            // slower code uses CMOV/select:
            //if (j.n != null) {
            //    s = if (i.v >= j.n.?.v) j else i;
            if (j.?.n != null and i.?.v >= j.?.n.?.v) {
                @branchHint(.likely);
                s = j.?;
            } else if (j.?.n != null) {
                @branchHint(.unlikely);
                s = i.?;
            } else if (j.?.v < h.?.v) {
                @branchHint(.unlikely);
                s = i.?;
            } else {
                @branchHint(.cold);
                self.reset();
                return [_]?*Node{null} ** 3;
            }

            const t = s.n.?;
            s.n = t.n;
            t.n = h.?;
            const iOut = if (t.v < h.?.v) t else i.?;

            return [_]?*Node{ t, iOut, iOut.n.? };
        }

        pub inline fn moveNext(self: *Self) bool {
            self.h, self.i, self.j = self.nextFrom(self.h, self.i, self.j);
            return self.h != null;
        }

        pub inline fn nextSeq(self: *Self) ?@Vector(16, u8) {
            if (self.h == null) {
                self.initNodes();
            } else if (!moveNext(self)) {
                return null;
            }

            return readLetters(self, self.h, self.input);
        }

        inline fn readLetters(self: *Self, head: ?*Node, input: @Vector(16, u8)) ?@Vector(16, u8) {
            var h = head;
            if (h == null) return null;

            const canMatchQuad = comptime std.meta.hasMethod(Visitor, "matchQuad");
            const Nquadded = 4 * (N / 4);

            var quads: @Vector(4, u32) = @bitCast(input);
            var x: u32 = 0;

            const canSkip = comptime std.meta.hasMethod(Visitor, "skipBeforeReading");
            if (canSkip and self.visitor.skipBeforeReading())
               return @bitCast(quads);
            
            inline for (0..Nquadded) |j| {
                const i = Nquadded - j - 1;
                x += h.?.r;
                h = h.?.n;

                // accumulate up to 4 bytes in x
                if ((i % 4) != 0) {
                    x = x << 8;
                    continue;
                }

                const q = 3 - j / 4;

                if (canMatchQuad) {
                    // bail out early if this quad fails
                    if (!self.visitor.matchQuad(x, q)) return null;
                }

                // update in quads
                quads[q] = x;
                if (i != 0) x = 0;
            }

            if (N == Nquadded)
                return @bitCast(quads);

            var letters: @Vector(16, u8) = @bitCast(quads);

            inline for (Nquadded..N) |j| {
                letters[MaxN - j - 1] = @intCast(h.?.r);
                h = h.?.n;
            }
            return letters;
        }

        pub fn nextMatch(self: *Self) ?@Vector(16, u8) {
            if (comptime !std.meta.hasMethod(Visitor, "match")) {
                @compileError("Visitor.match is not defined.");
            }

            // stuff that's faster if it's in a register
            const input: @Vector(16, u8) = self.input;

            if (self.h == null) {
                @call(.never_inline, initNodes, .{self});
                if (self.readLetters(self.h, input)) |puz| {
                    if (self.visitor.match(puz)) {
                        return puz;
                    }
                }
            }

            var h = self.h;
            var i = self.i;
            var j = self.j;

            while (true) {
                h, i, j = self.nextFrom(h, i, j);
                if (h == null) return null;

                if (self.readLetters(h, input)) |puz| {
                    if (self.visitor.match(puz)) {
                        self.h = h;
                        self.i = i;
                        self.j = j;
                        return puz;
                    }
                }
            }
        }
    };
}

var testCount: usize = 0;

fn CounterVisitorT() type {
    return struct {
        const Self = @This();

        pub inline fn skipBeforeReading(self: *const Self) bool {
            _ = self;
            testCount += 1;
            return true;
        }

        pub inline fn match(self: *const Self, x: @Vector(16, u8)) bool {
            _ = self;
            _ = x;
            return false;
        }
    };
}
const CounterVisitor = CounterVisitorT();

pub fn countMultiset(symbols: @Vector(16, u8), n: usize) usize {
    const visitor: CounterVisitor = .{};
    testCount = 0;

    inline for (1..(MaxN + 1)) |i| {
        if (i == n) {
            switch (i) {
                1 => {
                    return 1;
                },
                2 => {
                    return if (symbols[14] == symbols[15]) 1 else 2;
                },
                else => {
                    var mp: MultiPermuteWithVisitor(i, CounterVisitor) = .initVisitor(visitor);
                    mp.initLetters(symbols);
                    _ = mp.nextMatch();
                    return testCount;
                },
            }
        }
    }

    return 0;
}

pub const VisitorFn = *const fn (@Vector(16, u8)) void;

pub fn visitMultiset(symbols: @Vector(16, u8), n: usize, visit: VisitorFn) void {
    inline for (1..(MaxN + 1)) |i| {
        if (i == n) {
            switch (i) {
                1 => {
                    visit(symbols);
                },
                2 => {
                    visit(symbols);
                    if (symbols[14] != symbols[15]) {
                        var swapped = @Vector(16, u8){ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                        swapped = symbols;
                        swapped[14] = symbols[15];
                        swapped[15] = symbols[14];
                        visit(swapped);
                    }
                },
                else => {
                    var mp: MultiPermute(i) = .init();
                    mp.initLetters(symbols);
                    while (mp.nextSeq()) |seq| {
                        visit(seq);
                    }
                },
            }
        }
    }
}

test "counts" {
    var symbols = @Vector(16, u8){ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
    var count = countMultiset(symbols, 3);
    try testing.expectEqual(3, count);

    symbols = @Vector(16, u8){ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
    count = countMultiset(symbols, 4);
    try testing.expectEqual(4, count);

    symbols = @Vector(16, u8){ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 4 };
    count = countMultiset(symbols, 4);
    try testing.expectEqual(12, count);

    symbols = @Vector(16, u8){ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 1, 1 };
    count = countMultiset(symbols, 6);
    try testing.expectEqual(120, count);

    symbols = @Vector(16, u8){ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    count = countMultiset(symbols, 9);
    try testing.expectEqual(362880, count);
}
