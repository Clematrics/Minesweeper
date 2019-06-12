(* -------------------- Debug parameters -------------------- *)

let debug = false;;
let debug_visual = false;;

(* -------------------- Maths functions -------------------- *)

(* Returns the binomial coefficient i.e. p among n *)
let rec nCr n p =
	if p > n || p < 0 then
		0
	else if p = 0 || n = 1 then
		1
	else
		nCr (n - 1) (p - 1) * n / p
;;

(* Returns the binomial coefficient i.e. p among n with computation on 64 bits *)
let nCr64 n p =
	let rec aux n p =
		if Int64.equal p Int64.zero || Int64.equal n Int64.one then
			Int64.one
		else
			Int64.div (Int64.mul (aux (Int64.sub n Int64.one) (Int64.sub p Int64.one)) n) p
	in aux (Int64.of_int n) (Int64.of_int p)
;;

(* -------------------- Utils functions -------------------- *)

(* Prints a list and its content based on the format given, on the standard output *)
let print_list format l =
	print_string "[ ";
	List.iter (function i -> Printf.printf format i; print_string " ; ") l;
	print_string " ]";
	print_newline ()
;;

(* Prints an array and its content based on the func : 'a -> string,  on the standard output *)
let print_arr func l =
	print_string "[| ";
	Array.iter (function i -> Printf.printf "%s" (func i); print_string " ; ") l;
	print_string " |]";
	print_newline ()
;;

(* Makes 2^op operations that do nothing *)
let delay op = (* 27 or 28 is a good timing, to make a delay of a few seconds *)
	for i = 0 to 1 lsl op do
		()
	done
;;

(* -------------------- Array util functions -------------------- *)

(* Returns the dimension of a non-empty matrix *)
let dim arr =
	Array.length arr, Array.length arr.(0)
;;

(* Returns a matrix as an iterable list, line after line *)
let matrix_to_list matrix =
	let list_of_arr = Array.to_list matrix in
	List.concat (List.map (function l -> Array.to_list l) list_of_arr)
;;

(* Gives a list of neighbors of a cell (i,j) in the matrix of dimension (m, n) *)
let neighbors m n i j =
	let l = ref [] in
	if i > 0 then
		l := (i - 1, j):: !l;
	if j > 0 then
		l := (i, j - 1):: !l;
	if i > 0 && j > 0 then
		l := (i - 1, j - 1):: !l;
	if i > 0 && j < n - 1 then
		l := (i - 1, j + 1):: !l;
	if i < m - 1 then
		l := (i + 1, j):: !l;
	if j < n - 1 then
		l := (i, j + 1):: !l;
	if i < m - 1 && j > 0 then
		l := (i + 1, j - 1):: !l;
	if i < m - 1 && j < n - 1 then
		l := (i + 1, j + 1):: !l;
	!l
;;

(* -------------------- Generation of field -------------------- *)


(* Returns an array, which is indicator function of the seed-th configuration of count elements among total *)
let configuration total count seed =
	let rec aux l n p seed =
		if n = 0 then
			l
		else if p = 0 then
			aux (false::l) (n - 1) 0 seed
		else
			let proba = nCr (n - 1) (p - 1) in
			if seed < proba then
				aux (true::l) (n - 1) (p - 1) seed
			else
				aux (false::l) (n - 1) p (seed - proba)
	in
	if seed >= nCr total count then
		Printf.printf "Warning, seed too big for configuration generation";
	aux [] total count seed
;;

(* Same as above, but use the module Random to generate the field : good for generating fields *)
let configuration_opti total count seed =
	let conf = Array.make total false in
	let insert index =
		let i = ref index
		and j = ref 0 in
		while !i > -1 do
			if not conf.(!j) then
				decr i;
			incr j
		done;
		conf.(!j - 1) <- true
	in
	Random.init (seed);
	for i = 0 to count - 1 do
		insert (Random.int (total - i));
	done;
	Array.to_list conf
;;

(* Returns a field of dimension (m, n) with a number of landmines, with the (i,j) cell and its neighbors as safe cells, and with a seed *)
let field m n landmines i j seed =
	let neighbors = neighbors m n i j in
	let landmines_list = ref (configuration_opti (m * n - 1 - List.length neighbors) landmines seed)
	and field = Array.make_matrix m n false in
	for x = 0 to m - 1 do
		for y = 0 to n - 1 do
			if not (List.mem (x,y) neighbors) && (x,y) <> (i,j) then begin
				field.(x).(y) <- List.hd !landmines_list;
				landmines_list := List.tl !landmines_list
			end
		done
	done;
	field
;;

(* -------------------- Type descriptors -------------------- *)

type case = Hidden | Empty | Flag | Number of int | Landmine | Hypothesis of int * int;;
(* belongs to the grid which describes the current knowledge of the field *)
(* the hypothesis assumes an empty cell and guess the minimum and maximum number of landmines around it *)

type move = Dig of int * int | PutFlag of int * int | Explode of int * int;; (* describes a move *)

(* -------------------- Counting and functions -------------------- *)

(* Returns the number of landmines still hidden on the grid *)
let nb_mines_left grid total_mines =
	let m, n = dim grid in
	let coord_matrix = Array.init m (function i -> Array.init n (function j -> i,j)) in
	total_mines - List.fold_left (fun cnt -> function x,y -> cnt + (if grid.(x).(y) = Flag then 1 else 0)) 0 (matrix_to_list coord_matrix)
;;

(* Gives the number of landmines around a cell *)
let nb_landmines field i j =
	let m, n = Array.length field, Array.length field.(0) in
	let neighbors_ij = neighbors m n i j in
	List.fold_left	(fun cnt -> function x,y -> cnt + (if field.(x).(y) then 1 else 0)) 0 neighbors_ij
;;

(* Gives the number of hidden cells around a cell *)
let nb_hidden grid i j =
	let m, n = dim grid in
	let neighbors = neighbors m n i j in
	let hiddens = List.filter (function x,y -> grid.(x).(y) = Hidden) neighbors in
	List.length hiddens
;;

(* Gives the number of flags around a cell *)
let nb_flags grid i j =
	let m, n = dim grid in
	let neighbors = neighbors m n i j in
	let flags = List.filter (function x,y -> grid.(x).(y) = Flag) neighbors in
	List.length flags
;;

(* -------------------- Moves functions -------------------- *)

(* Expand the grid after a list of digs have been performed. Return the list of moves performed by the expansion *)
let expand grid field list =
	let rec aux dig_cells = function
		[] -> dig_cells
	|	(x,y)::t ->
			match grid.(x).(y) with
				Empty | Flag | Number _ -> aux dig_cells t
			|	_ ->
					let m, n = dim field in
					let cnt = nb_landmines field x y in
					if cnt = 0 then begin
						grid.(x).(y) <- Empty;
						let neighbors_xy = neighbors m n x y in
						let new_cells = List.filter
							(function x,y -> match grid.(x).(y) with
								Hidden -> true
							|	_ -> false
							) neighbors_xy in
						aux ((Dig (x,y))::dig_cells) (new_cells @ t)
					end
					else begin
						grid.(x).(y) <- Number cnt;
						aux ((Dig (x,y))::dig_cells) t
					end
	in aux [] list
;;

(* this is for real playing, should not return something *)
let apply_move grid field = function
	Dig (x, y) ->
		if grid.(x).(y) <> Hidden then begin
			if debug then
				Printf.printf "Warning, trying to dig in an inappropriate cell at %u, %u\n" x y;
			[]
		end
		else if field.(x).(y) then begin
			grid.(x).(y) <- Landmine;
			[Explode (x,y)]
		end
		else begin
			let landmines = nb_landmines field x y in
			if landmines = 0 then
				expand grid field [x, y]
			else begin
				grid.(x).(y) <- Number landmines;
				[Dig (x,y)]
			end
		end
|	PutFlag (x, y) ->
		if grid.(x).(y) = Flag then begin
			grid.(x).(y) <- Hidden;
			[PutFlag (x,y)]
		end
		else if grid.(x).(y) <> Hidden then begin
			if debug then 
				Printf.printf "Warning, trying to place a flag on an inappropriate cell at %u, %u\n" x y;
			[]
		end
		else begin
			grid.(x).(y) <- Flag;
			[PutFlag (x,y)]
		end
;;

let rec suppose_moves grid = function
	[] -> ()
|	Dig (x, y)::t ->
		if grid.(x).(y) <> Hidden then begin
			if debug then 
				Printf.printf "Warning, trying to dig in an inappropriate cell at %u, %u\n" x y
		end
		else
			grid.(x).(y) <- Hypothesis (1,8);
		suppose_moves grid t
|	PutFlag (x, y)::t ->
		if grid.(x).(y) <> Hidden then begin
			if debug then 
				Printf.printf "Warning, trying to place a flag on an inappropriate cell at %u, %u\n" x y
		end
		else
			grid.(x).(y) <- Flag;
		suppose_moves grid t
;;

let rec undo_moves grid = function
	[] -> ()
|	Dig (x, y)::t ->
		(match grid.(x).(y) with
			Number _ | Empty | Hypothesis _ -> grid.(x).(y) <- Hidden
			|	_ ->  if debug then Printf.printf "Warning, trying to undo digging on an inappropriate cell at %u, %u\n" x y
		);
		undo_moves grid t
|	PutFlag (x, y)::t ->
		if grid.(x).(y) <> Flag then begin
			if debug then 
				Printf.printf "Warning, trying to delete a flag on an inappropriate cell at %u, %u\n" x y
		end
		else
			grid.(x).(y) <- Hidden;
		undo_moves grid t
;;

(* -------------------- Grid generation -------------------- *)

(* Creates the grid based on the field and the initial dig cell *)
let init_grid field i j =
	let m, n = Array.length field, Array.length field.(0) in
	let grid = Array.make_matrix m n Hidden in
	grid.(i).(j) <- Empty;
	let neighbors_ij = neighbors m n i j in
	expand grid field neighbors_ij;
	grid
;;


(* -------------------- Display functions -------------------- *)

#load "graphics.cma";;

let cell_size = 60;;
let font_size = 20;;
let hidden_color = Graphics.rgb 128 128 128;;
let empty_color = Graphics.rgb 255 255 255;;
let flag_color = Graphics.rgb 200 20 20;;
let number_color = Graphics.rgb 196 196 196;;
let string_color = Graphics.white;;
let mine_color = Graphics.red;;

let display_flag x y =
	let flag = [| (2, 1) ; (2,5) ; (4,4) ; (2,3) |] in
	let poly = Array.map (function i,j -> x * cell_size + i * cell_size / 6, y * cell_size + j * cell_size / 6) flag in
	Graphics.set_color flag_color;
	Graphics.fill_poly poly
;;

let display_grid grid =
	let m, n = dim grid in
	for x = 0 to m - 1 do
		for y = 0 to n - 1 do
			let xlow = (x * cell_size)
			and xhigh = ((x + 1) * cell_size)
			and ylow = (y * cell_size)
			and yhigh = ((y + 1) * cell_size) in
			match grid.(x).(y) with
				Hidden ->
					Graphics.set_color hidden_color;
					Graphics.fill_rect xlow ylow cell_size cell_size
			|	Empty ->
					Graphics.set_color empty_color;
					Graphics.fill_rect xlow ylow cell_size cell_size
			|	Flag ->
					Graphics.set_color hidden_color;
					Graphics.fill_rect xlow ylow cell_size cell_size;
					display_flag x y
			|	Number n ->
					Graphics.set_color number_color;
					Graphics.fill_rect xlow ylow cell_size cell_size;
					Graphics.set_color string_color;
					Graphics.moveto ((xlow + xhigh) / 2 - font_size / 5 ) ((ylow + yhigh) / 2 - font_size / 3);
					Graphics.draw_string (Printf.sprintf "%u" n)
			|	Hypothesis _ ->
					Graphics.set_color number_color;
					Graphics.fill_rect xlow ylow cell_size cell_size;
					Graphics.set_color string_color;
					Graphics.moveto ((xlow + xhigh) / 2 - font_size / 5 ) ((ylow + yhigh) / 2 - font_size / 3);
					Graphics.draw_string (Printf.sprintf "?")
			|	Landmine ->
					Graphics.set_color mine_color;
					Graphics.fill_rect xlow ylow cell_size cell_size
		done
	done
;;

(* -------------------- IA functions -------------------- *);;

exception Abort_backtracking of move list list;;

(* Returns a move : place a flag or dig a hole. It also makes a guess *)
let evangeline_ia grid total_mines =
	let m, n = dim grid in
	let proba = Array.init m (function _ -> Array.init n (function _ -> [| 0 ; 0 |]) ) in (* the first value represent the number of configurations in which there is a mine, the second one represent the total of configurations *)
	let coord_matrix = Array.init m (function i -> Array.init n (function j -> i,j)) in
	let to_check = List.filter (function x, y ->
			match grid.(x).(y) with
				Number n -> nb_flags grid x y <= n && nb_hidden grid x y > 0
			|	_ -> false
		) (matrix_to_list coord_matrix) in
	(* backtracking evaluates the probability of a landmine in each cell *)
	let rec backtracking moves mines_left = function
		[] ->	let still_hidden = List.filter (function x,y -> grid.(x).(y) = Hidden) (matrix_to_list coord_matrix) in
				let nb_cells_left = List.length still_hidden in
				let nb_config = nCr nb_cells_left mines_left in
				List.iter (function
							Dig (x,y) ->
								proba.(x).(y).(1) <- nb_config + proba.(x).(y).(1);
						|	PutFlag (x,y) ->
								proba.(x).(y).(0) <- nb_config + proba.(x).(y).(0);
								proba.(x).(y).(1) <- nb_config + proba.(x).(y).(1);
					) (List.flatten moves);
				List.iter (function x,y ->
						proba.(x).(y).(0) <- nb_config * mines_left / nb_cells_left + proba.(x).(y).(0);
						proba.(x).(y).(1) <- nb_config + proba.(x).(y).(1)
						(*
							the number of configurations in which each unseen cell has a mine is nCr (nb_cells_left - 1) (mines_left - 1),
							but to avoid recomputing a binomial coefficient, we use a small property
						*)
					) still_hidden
	|	(x, y)::t ->
				if debug then
					Printf.printf "Dealing with %u,%u\n" x y;
				let neighbors = neighbors m n x y in
				let hidden_neighbors = List.filter (function x,y -> grid.(x).(y) = Hidden) neighbors
				and Number mines_around = grid.(x).(y) in
				let hidden_cnt = List.length hidden_neighbors
				and hidden_mines = mines_around - (nb_flags grid x y) in
				if hidden_cnt < hidden_mines then begin
					if debug then 
						Printf.printf "Contradiction : there are not enough hidden cells around to satisfy the landmines number!\n" (* Contradiction *)
				end
				else if hidden_mines > mines_left then begin
					if debug then 
						Printf.printf "Contradiction : there are not enough landmines on the field compared to hidden mines around a cell!\n" (* Contradiction *)
				end
				else if hidden_mines < 0 then begin
					if debug then 
						Printf.printf "Contradiction : There are too much flags around compared to the number of mines around!\n" (* Contradiction *)
				end
				else if hidden_mines = 0 then begin
					let new_moves = List.map (function x,y -> Dig (x,y)) hidden_neighbors in
					suppose_moves grid new_moves;
					backtracking (new_moves::moves) mines_left t;
					undo_moves grid new_moves
				end
				else begin
					for seed = 0 to nCr hidden_cnt hidden_mines - 1 do
						let configuration = configuration hidden_cnt hidden_mines seed in
						let combined = List.combine hidden_neighbors configuration in
						let flags = List.map (function coords, _ -> coords) (List.filter (function _, b -> b) combined)
						and digs = List.map (function coords, _ -> coords) (List.filter (function _, b -> not b) combined) in
						let flag_moves = List.map (function x,y -> PutFlag (x,y)) flags
						and dig_moves = List.map (function x,y -> Dig (x,y)) digs in
						suppose_moves grid flag_moves;
						suppose_moves grid dig_moves;
						(*delay 23;*)
						if debug_visual then
							display_grid grid;
						backtracking (flag_moves::dig_moves::moves) (mines_left - hidden_mines) t;
						undo_moves grid dig_moves;
						undo_moves grid flag_moves
 					done;
					(* List.iter (function x, y ->
						if ( proba.(x).(y).(0) = 0 && proba.(x).(y).(1) <> 0 ) || ( proba.(x).(y).(0) <> 0 && proba.(x).(y).(0) = proba.(x).(y).(1) ) then
							raise (Abort_backtracking moves) (* Moves can already be determined *)
					) hidden_neighbors *)
				end
	in
	(* try
		backtracking [] (nb_mines_left grid total_mines) to_check;
		raise (Abort_backtracking []); (* All probabilities where computed *)
		[Dig (0,0)] (* those two lines are just here to trick Caml thinking the type returned by the try-with is a move list like the function return type *)
	with
		Abort_backtracking moves -> undo_moves grid (List.flatten moves); *)

	backtracking [] (nb_mines_left grid total_mines) to_check;
	if debug then
		Array.iter (function arr -> print_arr (function [| a ; b |] -> Printf.sprintf "%u / %u" a b) arr) proba;
	let to_dig = List.filter (function x, y -> grid.(x).(y) = Hidden && proba.(x).(y).(0) = 0 && proba.(x).(y).(1) <> 0) (matrix_to_list coord_matrix)
	and to_flag = List.filter (function x, y -> grid.(x).(y) = Hidden && proba.(x).(y).(0) <> 0 && proba.(x).(y).(0) = proba.(x).(y).(1)) (matrix_to_list coord_matrix) in
	let dig_moves = List.map (function x, y -> Dig (x, y)) to_dig
	and flag_moves = List.map (function x, y -> PutFlag (x, y)) to_flag in
	let moves =	dig_moves @ flag_moves in
	if moves = [] then
		let to_manage = List.filter (function x, y -> proba.(x).(y).(1) <> 0) (matrix_to_list coord_matrix) in
		let x, y = List.fold_left (function minx, miny -> function x, y ->
				if proba.(x).(y).(0) < proba.(minx).(miny).(0) then
					x, y
				else
					minx, miny
			) (List.hd to_manage) (List.tl to_manage) in
		Printf.printf "I bet on the cell (%u, %u) not to hide a mine (probability is %f)\n" x y (float_of_int proba.(x).(y).(0) /. float_of_int proba.(x).(y).(1)) ;
		[Dig (x,y)]
	else
		moves
	(* TODO prise de decision du/des coups à faire : en faire plusieurs pour ce qui est sûr afin de gagner des appels au backtracking *)
;;

(* -------------------- Party managment -------------------- *)

let rec wait_move m n =
	let c = Graphics.read_key () in
	let x, y = Graphics.mouse_pos () in
	let i, j = x / cell_size, y / cell_size in
	if i < 0 || j < 0 || i >= m || j >= n then
		wait_move m n
	else if c = 'd' then
		Dig (i,j)
	else if c = 'f' then
		PutFlag (i,j)
	else
		wait_move m n
;;

let party m n landmines i j seed =
	Graphics.open_graph (Printf.sprintf "%ux%u+500+500" (cell_size * m) (cell_size * n));
	Graphics.set_text_size font_size;
	let field = field m n landmines i j seed in
	let grid = init_grid field i j in
	let playing = ref true in
	while !playing do
		display_grid grid;
		let move = wait_move m n in
		let res = apply_move grid field move in
		match res with
			Explode (x,y)::_ ->
				Printf.printf "Game over!\n";
				playing := false
		|	_ ->
				let coord_matrix = Array.init m (function i -> Array.init n (function j -> i,j)) in
				let good = List.for_all (function x, y -> if grid.(x).(y) = Hidden || grid.(x).(y) = Flag then field.(x).(y) else true) (matrix_to_list coord_matrix) in
				if good then begin
					Printf.printf "You won!\n";
					playing := false
				end
	done;
	display_grid grid;
;;

exception Win;;
exception GameOver;;

let party_ia m n landmines i j seed =
	Graphics.open_graph (Printf.sprintf "%ux%u+500+500" (cell_size * m) (cell_size * n));
	Graphics.set_text_size font_size;
	let field = field m n landmines i j seed in
	let grid = init_grid field i j in
	let coord_matrix = Array.init m (function i -> Array.init n (function j -> i,j)) in
	try
		while true do
			display_grid grid;
			let good = List.for_all (function x, y -> if grid.(x).(y) = Hidden || grid.(x).(y) = Flag then field.(x).(y) else true) (matrix_to_list coord_matrix) in
			if good then
				raise Win;
			let moves = evangeline_ia grid landmines in
			List.iter (function move ->
				let res = apply_move grid field move in
				match res with
					Explode (_)::_ -> raise GameOver
				|	_ -> ()
			) moves;
		done;
		true (* line to trick OCaml *)
	with
		Win ->
			Printf.printf "You win!\n";
			display_grid grid;
			true
	|	GameOver ->
			Printf.printf "Game over!\n";
			display_grid grid;
			false
;;

party 10 10 20 0 0 73788;;

party_ia 20 20 50 0 0 73788;;

let cnt = ref 0 in
for i = 0 to 20 do
	let won = party_ia 10 16 34 5 8 i in
	if won then incr cnt;
done;
Printf.printf "%u" !cnt
;;