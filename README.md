# MultiPermute.zig

Efficient multiset permutations,
adapted from https://github.com/ekg/multipermute
and optimized for speed and WebAssembly compatibility.

Does not work with sequences longer than 1024.

## Overview

This package contains a method to generate all permutations of a multiset. The method is described in ["Loopless Generation of Multiset Permutations using a Constant Number of Variables by Prefix Shifts."](https://www.researchgate.net/publication/220779544_Loopless_Generation_of_Multiset_Permutations_using_a_Constant_Number_of_Variables_by_Prefix_Shifts) Aaron Williams, 2009. To quote [Aaron's website](http://webhome.cs.uvic.ca/~haron/):

> The permutations of any multiset are generated by this simple rule:

> - The first symbol a is shifted to the right until it passes over consecutive symbols b c with b < c. 
> - If a > b, then a is inserted after b; otherwise, if a <= b, then a is inserted after c. 
> - (If there is no such b c then a is shifted until it passes over the rightmost symbol.)

> This result leads to the first O(1)-time / O(1)-additional variable algorithm for generating the permutations of a multiset (see publication in SODA 2009)."

## Usage

``` bash
% zig build -Doptimize=ReleaseFast
```

```bash
% zig-out/bin/mp_list 0 1 2 2
2210
0221
2021
2201
1220
2120
0212
2012
1202
0122
1022
2102
```

```bash
% zig-out/bin/mp_count 0 1 2 2
12
```

### Basic Zig API

```
zig fetch --save 'https://github.com/nicolaitanes/MultiPermute.zig/archive/refs/tags/v1.0.0.tar.gz'
```

```zig
const MultiPermute = @import("MultiPermute");

pub fn main() {
    const N = 12;
    var symbols: [_]u8 = .{65,65,66,66,67,67,67,68,68,69,70,71};
    
    const count = MultiPermute.countMultiset(symbols[0..N]);
    
    MultiPermute.visitMultiset(symbols[0..N], &myVisitor);
}

fn myVisitor(seq: []const u8) void {
    // do something with seq...
}
```

### Full Zig API

#### Types

- `MultiPermute()`: generates permutations of up to 1024 symbols
- `MultiPermuteWithVisitor(Visitor)`: uses a visitor to find or count matching permutations
- `MultiPermuteSmall(N)`: keeps the whole sequence in a register (avoids memory accesses), when sequence length `N` is known at compile time and between 3 and 16 (inclusive),
- `MultiPermuteSmallWithVisitor(N, Visitor)`

#### Module import

```
zig fetch --save 'https://github.com/nicolaitanes/MultiPermute.zig/archive/refs/tags/v1.0.0.tar.gz'
```

zig.build:

```
    const package = b.dependency("MultiPermute", .{
        .target = target,
        .optimize = optimize,
    });
    const module = package.module("MultiPermute")

    myexe.root_module.addImport("MultiPermute", module);
```

#### Iterating through permutations

```zig
const MultiPermute = @import("MultiPermute");

pub fn main() {
    const N = 12;
    var symbols: [_]u8 = .{65,65,66,66,67,67,67,68,68,69,70,71};

    var mp12: MultiPermute() = .init();
    mp12.initLetters(symbols[0..N]);

    while (mp12.nextSeq()) |seq| {
        // do something with seq...
    }
}
```

#### Iterating through permutations (3 <= N <= 16)

```zig
const MultiPermute = @import("MultiPermute");

pub fn main() {
    const N = 12;
    var symbols = @Vector(16, u8){0,0,0,0,65,65,66,66,67,67,67,68,68,69,70,71};

    var mp12: MultiPermuteSmall(N) = .init();
    mp12.initLetters(symbols);

    while (mp12.nextSeq()) |seq| {
        // do something with seq...
    }
}
```

#### Search with custom Visitor (3 <= N <= 16)

Symbol order is kept in a linked list; the algorithm can spend a significant amount of time traversing the nodes to read each permutation.

If some or all output sequences are not actually needed (e.g. when counting permutations), use `MultiPermute[Small]WithVisitor` to short-circuit list traversal:

```zig
const MultiPermute = @import("MultiPermute");

fn VisitorT() type {
    return struct {
        const Self = @This();

        // *optional* return true to skip all list traversal (e.g. for counting).
        pub inline fn skipBeforeReading(self: *const Self) bool {
            return shouldBeSkipped(x, i);
        }

        // *optional* allows inspection of each completed set of 4 characters in `x` (starting from i=3, down to i=0);
        // return false to skip to the next permutation.
        pub inline fn matchQuad(self: *const Self, x: u32, comptime i: usize) bool {
            return shouldBeSkipped(x, i);
        }

        // *required* if calling `MultiPermuteWithVisitor.nextMatch`;
        // return true if the completed permutation in `x` should be returned from nextMatch
        // (pausing iteration)
        pub inline fn match(self: *const Self, x: @Vector(16, u8)) bool {
            return isMatching(x);
        }
    };
}
const Visitor = VisitorT();

pub fn main() {
    const N = 12;
    var symbols = @Vector(16, u8){0,0,0,0,65,65,66,66,67,67,67,68,68,69,70,71};

    var mp12: MultiPermuteSmallWithVisitor(N, Visitor) = .initVisitor(.{});
    mp12.initLetters(symbols);

    while (mp12.nextMatch()) |seq| {
        // do something with seq...
    }
}
```
