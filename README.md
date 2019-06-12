# Minesweeper

A group of projects about the minesweeper

# Implementation in OCaml

A minesweeper game and AI implemented in OCaml.

You can play yourself by calling the function `party m n landmines i j seed`, where :
  * `m` is the width of the field
  * `n` is the height of the field
  * `landmines` is the number of landmines you want
  * `i, j` is the coordinates of the initial move (so you don't trigger a mine right away)
  * `seed` is the seed of the field

To dig a cell, place your mouse above it and press `d`.
To place or take off a flag, place your mouse above the cell and press `f`.


You can also watch an AI solve the game by calling the function `party_ai m n landmines i j seed`. The parameters are the same as above.
Because OCaml cannot handle too big numbers, do not use too large fields (30x30 is a maximum), or the AI will take forever to solve the game.
The AI uses backtracing to compute probabilities of the presence of a mine on every hidden cell. It detects contradictions and make the best moves according to the probabilities. However, the AI does not take yet into account information it can obtain from a cell, despite the probability is not minimal, so there is room for improvement.

![OCaml screenshot](https://github.com/Clematrics/Minesweeper/blob/master/ocaml%20screenshot.png)