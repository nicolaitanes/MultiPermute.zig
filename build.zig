const std = @import("std");
const Target = @import("std").Target;
const Feature = @import("std").Target.Cpu.Feature;

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const mp_count = b.addExecutable(.{ .name = "mp_count", .root_source_file = b.path("src/mp_count.zig"), .target = target, .optimize = optimize, .single_threaded = true, .strip = true });
    const mp_list = b.addExecutable(.{ .name = "mp_list", .root_source_file = b.path("src/mp_list.zig"), .target = target, .optimize = optimize, .single_threaded = true, .strip = true });

    b.installArtifact(mp_count);
    b.installArtifact(mp_list);
}
