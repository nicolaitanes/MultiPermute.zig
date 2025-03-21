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

pub const MaxN = 1024;

/// This module encodes functions to generate the permutations of a multiset
/// following this algorithm:
///
/// Algorithm 1
/// Visits the permutations of multiset E. The permutations are stored
/// in a singly-linked list pointed to by head pointer h. Each node in the linked
/// list has a value field v and a next field n. The init(E) call creates a
/// singly-linked list storing the elements of E in non-increasing order with h, i,
/// and j pointing to its first, second-last, and last nodes, respectively. The
/// null pointer is given by φ. Note: If E is empty, then init(E) should exit.
/// Also, if E contains only one element, then init(E) does not need to provide a
/// value for i.
///
/// A linked list must be used instead of an array because, as long as we keep a
/// reference to the nodes at the removal and insertion points, arbitrary shifts
/// (rotations of sub-lists) can be performed in constant time; this is what h, i,
/// and j are for.
///
/// [h, i, j] ← init(E)
/// visit(h)
/// while j.n ≠ φ or j.v < h.v do
///     if j.n ≠ φ and i.v ≥ j.n.v then
///         s←j
///     else
///         s←i
///     end if
///     t←s.n
///     s.n ← t.n
///     t.n ← h
///     if t.v < h.v then
///         i←t
///     end if
///     j←i.n
///     h←t
///     visit(h)
/// end while
///
/// Adapted from https://github.com/ekg/multipermute
/// and optimized for fast WebAssembly, with a limit of 3-16 elements in the multiset.
pub fn MultiPermute() type {
    return MultiPermuteWithVisitor(struct {});
}

pub const MaxNsmall = 16;

// Specialized to keep results in a register, avoiding ram.
pub fn MultiPermuteSmall(comptime N: u32) type {
    return MultiPermuteSmallWithVisitor(N, struct {});
}

fn ms2permutation(symbols: []const u8, permu: []u8, remap: []u8, N: u16) void {
    // ms2multiplicities
    const capacity = symbols.len;
    var counters = [_]u8{0} ** 256;
    for (0..N) |i| {
        const e = symbols[capacity - i - 1];
        counters[e] += 1;
    }
    var mults = [_]u8{0} ** MaxN;
    var nMult: usize = 0;
    for (0..N) |i| {
        const e = symbols[capacity - i - 1];
        const n = counters[e];
        if (n == 0) continue;
        remap[nMult] = e;
        mults[nMult] = n;
        counters[e] = 0;
        nMult += 1;
    }
    // multiplicity2sorted:
    var j: usize = 0;
    for (0..N) |k| {
        var m = mults[k];
        while (m > 0) : (m -= 1) {
            permu[j] = @intCast(k);
            j += 1;
        }
    }
}

const MultiPermuteNode = struct { v: u32, r: u32, n: ?*MultiPermuteNode };

fn initMultiPermuteList(nodes: []MultiPermuteNode, permu: []u8, remap: []u8) [3]?*MultiPermuteNode {
    const N = nodes.len;
    for (0..N) |i| {
        nodes[i].v = permu[N - i - 1];
        const v: u32 = @bitCast(nodes[i].v);
        nodes[i].r = remap[v];
        nodes[i].n = if (i < N - 1) &nodes[i + 1] else null;
    }

    return [_]?*MultiPermuteNode{ &nodes[0], &nodes[N - 2], &nodes[N - 1] };
}

fn readMultiPermuteList(head: ?*MultiPermuteNode, symbols: []u8) bool {
    var h = head;
    if (h == null) return false;

    for (0..symbols.len) |j| {
        const i = symbols.len - j - 1;
        symbols[i] = @intCast(h.?.r);
        h = h.?.n;
    }

    return true;
}

/// Rewires the linked list [h, ..., i, j] for the next permutation, and returns the next [h, i, j]
// ([head, penultimate, ultimate]).
inline fn nextFrom(h: ?*MultiPermuteNode, i: ?*MultiPermuteNode, j: ?*MultiPermuteNode) [3]?*MultiPermuteNode {
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
        return [_]?*MultiPermuteNode{null} ** 3;
    }

    const t = s.n.?;
    s.n = t.n;
    t.n = h.?;
    const iOut = if (t.v < h.?.v) t else i.?;

    return [_]?*MultiPermuteNode{ t, iOut, iOut.n.? };
}

/// For applications where reading each permutation (traversing a linked list)
/// is a significant fraction of runtime, provide a Visitor that can matchQuad(fourBytes, index in 0..4) --
/// if not matched, list traversal is abandoned and we move on to the next permutation.
/// Visitor: { .skipBeforeReading() bool, .matchQuad(x, i) bool, .match(x) bool }
pub fn MultiPermuteWithVisitor(comptime Visitor: type) type {
    return struct {
        N: u16 = 0,
        visitor: Visitor,
        nodes: [MaxN]Node,
        permu: [MaxN]u8,
        remap: [MaxN]u8,
        h: ?*Node = null,
        i: ?*Node = null,
        j: ?*Node = null,
        input: [MaxN]u8,

        const Self = @This();
        const Node = MultiPermuteNode;

        pub inline fn initVisitor(visitor: Visitor) Self {
            return Self{ .visitor = visitor, .nodes = [_]Node{Node{ .v = 0, .r = 0, .n = null }} ** MaxN, .permu = [_]u8{0} ** MaxN, .remap = [_]u8{0} ** MaxN, .input = [_]u8{0} ** MaxN };
        }

        pub inline fn init() Self {
            return initVisitor(.{});
        }

        pub inline fn reset(self: *Self) void {
            self.h = null;
        }

        /// Inputs symbols for permutation.
        pub fn initSymbols(self: *Self, symbols: []const u8) void {
            self.N = @min(symbols.len, MaxN);
            for (0..self.N) |i|
                self.input[i] = symbols[i];

            ms2permutation(symbols, &self.permu, &self.remap, self.N);

            self.reset();
        }

        /// Initializes the linked list with input symbols and initial traversal order.
        pub fn initNodes(self: *Self) void {
            if (self.permu[self.N - 1] == 0)
                return self.initSymbols(self.input[0..self.N]);

            self.h, self.i, self.j = initMultiPermuteList(self.nodes[0..self.N], self.permu[0..self.N], self.remap[0..self.N]);
        }

        /// Moves to the next permutation.
        pub inline fn moveNext(self: *Self) bool {
            self.h, self.i, self.j = nextFrom(self.h, self.i, self.j);
            if (self.h == null) {
                self.reset();
                return false;
            }
            return true;
        }

        /// Moves to the first or next permutation (if any) and reads it out.
        pub inline fn nextSeq(self: *Self) ?[]u8 {
            if (self.h == null) {
                self.initNodes();
            } else if (!moveNext(self)) {
                return null;
            }

            return readSymbols(self, self.h);
        }

        /// Reads out the current permutation, or null if there are no more.
        inline fn readSymbols(self: *Self, head: ?*Node) ?[]u8 {
            var h = head;
            if (h == null) return null;

            const canMatchQuad = comptime std.meta.hasMethod(Visitor, "matchQuad");
            const Nquadded = 4 * (self.N / 4);
            var x: u32 = 0;

            const canSkip = comptime std.meta.hasMethod(Visitor, "skipBeforeReading");
            if (canSkip and self.visitor.skipBeforeReading())
                return self.input[0..self.N];

            for (0..self.N) |j| {
                const r: u32 = h.?.r;
                self.input[self.N - j - 1] = @intCast(r);
                if (j < (self.N - 1))
                    h = h.?.n;

                if (canMatchQuad and j < Nquadded) {
                    x += r;
                    const i = Nquadded - j - 1;
                    if ((i % 4) != 0) {
                        x = x << 8;
                    } else if (!self.visitor.matchQuad(x, 3 - j / 4)) {
                        return null;
                    } else if (i != 0) {
                        x = 0;
                    }
                }
            }

            return self.input[0..self.N];
        }

        /// Iterates through permutations, returning the first one for which `Visitor.match(seq)` is true.
        pub fn nextMatch(self: *Self) ?[]u8 {
            if (comptime !std.meta.hasMethod(Visitor, "match")) {
                @compileError("Visitor.match is not defined.");
            }

            if (self.h == null) {
                @call(.never_inline, initNodes, .{self});
                if (self.readSymbols(self.h)) |puz| {
                    if (self.visitor.match(puz)) {
                        return puz;
                    }
                }
            }

            var h = self.h;
            var i = self.i;
            var j = self.j;

            while (true) {
                h, i, j = nextFrom(h, i, j);
                if (h == null) {
                    self.reset();
                    return null;
                }

                if (self.readSymbols(h)) |puz| {
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

/// For applications where reading each permutation (traversing a linked list)
/// is a significant fraction of runtime, provide a Visitor that can matchQuad(fourBytes, index in 0..4) --
/// if not matched, list traversal is abandoned and we move on to the next permutation.
/// Visitor: { .skipBeforeReading() bool, .matchQuad(x, i) bool, .match(x) bool }
pub fn MultiPermuteSmallWithVisitor(comptime N: u32, comptime Visitor: type) type {
    if (N > MaxNsmall) @compileError("Max 16 symbols.");
    return struct {
        visitor: Visitor,
        nodes: [N]Node,
        permu: [N]u8,
        remap: [N]u8,
        N: u32 = N,
        h: ?*Node = null,
        i: ?*Node = null,
        j: ?*Node = null,
        input: @Vector(16, u8) = @splat(0),

        const Self = @This();
        const Node = MultiPermuteNode;

        pub inline fn initVisitor(visitor: Visitor) Self {
            return Self{ .visitor = visitor, .nodes = [_]Node{Node{ .v = 0, .r = 0, .n = null }} ** N, .permu = [_]u8{0} ** N, .remap = [_]u8{0} ** N };
        }

        pub inline fn init() Self {
            return initVisitor(.{});
        }

        pub inline fn reset(self: *Self) void {
            self.h = null;
        }

        /// Inputs symbols for permutation. The first (N-16) symbols will be untouched in the output;
        /// the remaining N symbols will be permuted.
        pub fn initSymbols(self: *Self, symbols: @Vector(16, u8)) void {
            self.input = symbols;

            var symbolArr = [_]u8{0} ** 16;
            symbolArr[0..].* = symbols;
            ms2permutation(&symbolArr, &self.permu, &self.remap, N);

            self.reset();
        }

        /// Initializes the linked list with input symbols and initial traversal order.
        pub fn initNodes(self: *Self) void {
            if (self.permu[N - 1] == 0)
                return self.initSymbols(self.input);

            self.h, self.i, self.j = initMultiPermuteList(&self.nodes, &self.permu, &self.remap);
        }

        /// Moves to the next permutation.
        pub inline fn moveNext(self: *Self) bool {
            self.h, self.i, self.j = nextFrom(self.h, self.i, self.j);
            if (self.h == null) {
                self.reset();
                return false;
            }
            return true;
        }

        /// Moves to the first or next permutation (if any) and reads it out.
        pub inline fn nextSeq(self: *Self) ?@Vector(16, u8) {
            if (self.h == null) {
                self.initNodes();
            } else if (!moveNext(self)) {
                return null;
            }

            return readSymbols(self, self.h, self.input);
        }

        /// Reads out the current permutation, or null if there are no more.
        inline fn readSymbols(self: *Self, head: ?*Node, input: @Vector(16, u8)) ?@Vector(16, u8) {
            var h = head;
            if (h == null) return null;

            const canMatchQuad = comptime std.meta.hasMethod(Visitor, "matchQuad");
            const Nquadded = 4 * (N / 4);

            var symbols: @Vector(16, u8) = @bitCast(input);
            var x: u32 = 0;

            const canSkip = comptime std.meta.hasMethod(Visitor, "skipBeforeReading");
            if (canSkip and self.visitor.skipBeforeReading())
                return input;

            inline for (0..N) |j| {
                const r: u32 = h.?.r;
                symbols[MaxNsmall - j - 1] = @intCast(r);
                if (j < (N - 1))
                    h = h.?.n;

                // accumulate up to 4 bytes in x
                if (canMatchQuad and j < Nquadded) {
                    x += r;
                    const i = Nquadded - j - 1;
                    if ((i % 4) != 0) {
                        x = x << 8;
                    } else if (!self.visitor.matchQuad(x, 3 - j / 4)) {
                        return null;
                    } else if (i != 0) {
                        x = 0;
                    }
                }
            }

            return symbols;
        }

        /// Iterates through permutations, returning the first one for which `Visitor.match(seq)` is true.
        pub fn nextMatch(self: *Self) ?@Vector(16, u8) {
            if (comptime !std.meta.hasMethod(Visitor, "match")) {
                @compileError("Visitor.match is not defined.");
            }

            // stuff that's faster if it's in a register
            const input: @Vector(16, u8) = self.input;

            if (self.h == null) {
                @call(.never_inline, initNodes, .{self});
                if (self.readSymbols(self.h, input)) |puz| {
                    if (self.visitor.match(puz)) {
                        return puz;
                    }
                }
            }

            var h = self.h;
            var i = self.i;
            var j = self.j;

            while (true) {
                h, i, j = nextFrom(h, i, j);
                if (h == null) {
                    self.reset();
                    return null;
                }

                if (self.readSymbols(h, input)) |puz| {
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

        pub inline fn match(self: *const Self, x: []u8) bool {
            _ = self;
            _ = x;
            return false;
        }
    };
}
const CounterVisitor = CounterVisitorT();

// convenience wrapper for counting permutations

/// Returns the number of distinct multisets that can be made from `symbols`.
pub fn countMultiset(symbols: []const u8) usize {
    switch (symbols.len) {
        0...1 => {
            return 1;
        },
        2 => {
            return if (symbols[0] == symbols[1]) 1 else 2;
        },
        else => {
            const visitor: CounterVisitor = .{};
            testCount = 0;

            var mpLarge: MultiPermuteWithVisitor(CounterVisitor) = .initVisitor(visitor);
            mpLarge.initSymbols(symbols);
            _ = mpLarge.nextMatch();

            return testCount;
        },
    }
}

// convenience wrapper for visiting each permutation

pub const VisitorFn = *const fn ([]const u8) void;

/// Calls `visit` on each distinct multiset that can be made from the right-most `n` items in `symbols`.
pub fn visitMultiset(symbols: []const u8, visit: VisitorFn) void {
    var symbolIO = [_]u8{0} ** MaxN;
    switch (symbols.len) {
        0...1 => {
            visit(symbols);
        },
        2 => {
            visit(symbols);
            if (symbols[0] != symbols[1]) {
                symbolIO[0] = symbols[1];
                symbolIO[1] = symbols[0];
                visit(symbolIO[0..]);
            }
        },
        else => {
            var mpLarge: MultiPermute() = .init();
            mpLarge.initSymbols(symbols);
            while (mpLarge.nextSeq()) |seq| {
                visit(seq);
            }
        },
    }
}
