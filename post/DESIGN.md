# PostScript — G-code Post Processor Language

A turing-complete, Python-style scripting language for generating G-code.
Event-driven: the post engine walks the toolpath and calls your handlers.

## Why not a template language

Every "simple" post template language eventually needs loops, conditionals,
state tracking, arithmetic, and string formatting. Then it becomes a bad
programming language. Mastercam's post processor is the cautionary tale —
pattern matching that grew into an unmaintainable mess of IF/THEN/ELSE
blocks, buffer variables, and post flags.

Start with a real language. Keep it small and readable.

## Design Principles

1. **Python-style syntax** — indent-scoped blocks, no braces, no semicolons
2. **Event-driven** — you define handlers, the engine calls them
3. **F-string formatting** — `"G01 X{move.x:.4}"` for output
4. **Full control flow** — if/elif/else, for loops, while loops, functions
5. **CAM-aware standard library** — knows about tools, coordinates, feeds, arcs
6. **Tree-sitter grammar** — syntax highlighting, auto-indent, error squiggles
7. **LLM-writable** — "set up a post for my Haas VF-4" → working post

## Syntax Example

```
# Haas VF-4 — 3-axis vertical mill
machine "Haas VF-4"
controller "Haas NGC"

def on_program_start(program):
    output "%"
    output "O{program.number:04d} ({program.name})"
    output "G20 G90 G40 G80"    # inch, absolute, cancel comp, cancel canned

def on_tool_change(tool, next_op):
    output "T{tool.number} M06 ({tool.description})"
    output "G43 H{tool.number} Z{tool.gauge_length:.4}"
    output "S{next_op.speed} M03"

    # Shop-specific: turn coolant on 1" above feed plane during approach,
    # not at tool change. See on_rapid_move.

def on_rapid_move(move):
    # Coolant-on logic: activate 1" above the feed plane on Z approach.
    # This prevents coolant from running while the tool is way up high
    # (keeps the window clean, saves coolant on long tools).
    if move.axis == 'Z' and move.direction == 'down':
        feed_plane = current_op.feed_plane
        coolant_on_at = feed_plane + 1.0
        if move.start.z > coolant_on_at and move.end.z <= coolant_on_at:
            output "G00 Z{coolant_on_at:.4}"
            if current_op.tool.coolant == "through":
                output "M88"
            else:
                output "M08"
            output "G00 Z{move.end.z:.4}"
            return

    output "G00 {move.format()}"

def on_linear_move(move):
    output "G01 X{move.x:.4} Y{move.y:.4} Z{move.z:.4} F{move.feed:.1}"

def on_arc_cw(arc):
    if arc.is_full_circle:
        # Haas can't do full circles in one block — split at 180°
        mid = arc.midpoint()
        output "G02 X{mid.x:.4} Y{mid.y:.4} R{arc.radius:.4} F{arc.feed:.1}"
        output "G02 X{arc.x:.4} Y{arc.y:.4} R{arc.radius:.4}"
    else:
        output "G02 X{arc.x:.4} Y{arc.y:.4} R{arc.radius:.4} F{arc.feed:.1}"

def on_arc_ccw(arc):
    if arc.is_full_circle:
        mid = arc.midpoint()
        output "G03 X{mid.x:.4} Y{mid.y:.4} R{arc.radius:.4} F{arc.feed:.1}"
        output "G03 X{arc.x:.4} Y{arc.y:.4} R{arc.radius:.4}"
    else:
        output "G03 X{arc.x:.4} Y{arc.y:.4} R{arc.radius:.4} F{arc.feed:.1}"

def on_coolant_off():
    output "M09"

def on_spindle_stop():
    output "M05"

def on_program_end(program):
    output "M05"
    output "M09"
    output "G91 G28 Z0"
    output "G28 X0 Y0"
    output "M30"
    output "%"
```

## The Coolant Example

"Turn the coolant on right at 1 inch from the feed plane."

This one sentence is the litmus test. It requires:
- Knowing the current Z position
- Knowing the operation's feed plane Z
- Calculating when you're 1" above it during the Z approach
- Splitting a rapid move into two segments
- Inserting M08/M88 between them

In Mastercam's post, this is a rat's nest of buffer variables and post flags.
In PostScript, it's 8 lines of readable code in `on_rapid_move`. A machinist
can read it, understand it, and change "1 inch" to "0.5 inch" without calling
a post consultant.

## Event Model

The post engine walks the toolpath and fires events:

| Event | When | Data |
|-------|------|------|
| `on_program_start` | Beginning of program | program name, number |
| `on_program_end` | End of program | program name |
| `on_tool_change` | Tool change | tool spec, next operation |
| `on_operation_start` | Start of an operation | operation params |
| `on_operation_end` | End of an operation | operation params |
| `on_rapid_move` | G00 rapid positioning | start, end, axis |
| `on_linear_move` | G01 feed move | x, y, z, feed |
| `on_arc_cw` | G02 clockwise arc | endpoint, radius/center, feed |
| `on_arc_ccw` | G03 counterclockwise arc | endpoint, radius/center, feed |
| `on_drill_cycle` | G81-G89 canned cycle | cycle type, params |
| `on_coolant_on` | Coolant activation | type (flood/mist/through) |
| `on_coolant_off` | Coolant deactivation | |
| `on_spindle_start` | Spindle on | RPM, direction |
| `on_spindle_stop` | Spindle off | |
| `on_dwell` | G04 dwell | duration |
| `on_comment` | Insert comment | text |

If you don't define a handler, the engine uses a sensible default.
Override only what your machine needs.

## Built-in Objects

```
current_op          # current operation being posted
current_op.tool     # tool for this operation
current_op.feed_plane   # Z of the feed plane (retract height)
current_op.speed    # RPM
current_op.feed     # IPM

move.x, move.y, move.z   # endpoint coordinates
move.start          # start position (Point3)
move.end            # end position (Point3)
move.axis           # dominant axis ('X', 'Y', 'Z')
move.direction      # 'up' or 'down' (for Z moves)
move.feed           # feed rate (IPM)
move.format()       # auto-format changed axes only (modal suppression)

arc.center          # arc center point
arc.radius          # arc radius
arc.start_angle     # start angle
arc.end_angle       # end angle
arc.is_full_circle  # true if 360°
arc.midpoint()      # point at 180° (for full-circle splitting)

tool.number         # tool number
tool.diameter       # tool diameter
tool.description    # tool description string
tool.coolant        # "flood", "mist", "through", "off"
tool.gauge_length   # gauge length
tool.flutes         # number of flutes
```

## Editor Integration

- Tree-sitter grammar for syntax highlighting and structure
- Auto-indent on newline (Python-style block detection)
- Error squiggles from the interpreter's parse phase
- Autocomplete for event names, built-in objects, and properties
- Built into the ForgeCAM GUI as a tab in the workspace panel

## Anti-Patterns (What We're NOT Doing)

- **NOT a template language** — no `%VARIABLE%` substitution, no pattern matching
- **NOT Lua** — too general, no indent-scoping, the machinist shouldn't need
  to learn metatables
- **NOT embedded Python** — too heavy a runtime dependency for what we need,
  and we want the tree-sitter grammar to be ours
- **NOT Mastercam-style** — no buffer variables, no implicit state machine,
  no "post flags" that silently change behavior

## Status

Design phase. No code yet. Will live in `post/` crate.
