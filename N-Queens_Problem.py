import tkinter as tk
from tkinter import messagebox
import random
import colorsys

# -------------------------------
# Colors
# -------------------------------

def generate_palette(n):
    """Distinct-ish colors via HSV wheel."""
    out = []
    for k in range(n):
        h = k / max(1, n)
        s = 0.55
        v = 0.95
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        out.append('#{:02x}{:02x}{:02x}'.format(int(r*255), int(g*255), int(b*255)))
    return out

# -------------------------------
# Region generator: contiguous irregular blobs
# -------------------------------

def generate_blob_regions(n, max_restarts=200):
    """
    Generate N contiguous, irregular blobs that cover an n×n grid.
    Each region is a single connected component (4-neighborhood).
    If the growth gets stuck, we restart (up to max_restarts).
    """
    dirs = [(1,0), (-1,0), (0,1), (0,-1)]

    for _attempt in range(max_restarts):
        board = [[-1]*n for _ in range(n)]
        target_sizes = [n]*n  # total cells = n*n, so each region gets n cells
        ok = True

        # grow each region fully before moving to the next (keeps blobs compact)
        for rid in range(n):
            # pick a random unassigned seed
            unassigned = [(i,j) for i in range(n) for j in range(n) if board[i][j] == -1]
            if not unassigned:
                ok = False
                break
            sr, sc = random.choice(unassigned)
            board[sr][sc] = rid
            size = 1
            frontier = [(sr, sc)]

            while size < target_sizes[rid]:
                # collect all frontier-adjacent unassigned neighbors
                candidates = []
                for fr, fc in frontier:
                    for di, dj in dirs:
                        nr, nc = fr+di, fc+dj
                        if 0 <= nr < n and 0 <= nc < n and board[nr][nc] == -1:
                            candidates.append((nr, nc))
                if not candidates:
                    ok = False  # this region cannot grow to its target -> restart
                    break

                nr, nc = random.choice(candidates)
                board[nr][nc] = rid
                frontier.append((nr, nc))
                size += 1

            if not ok:
                break

        if ok:
            return board

    # If we failed to generate after many restarts, fall back to a simple split (still contiguous)
    # but in practice the loop above should succeed quickly for typical N (≤ 14).
    board = [[(i // max(1, n//n)) % n for j in range(n)] for i in range(n)]
    return board

# -------------------------------
# N-Queens with color constraint
# -------------------------------

def can_place_manual(board, row, col, n, region):
    # row
    for j in range(n):
        if board[row][j] == 1:
            return False
    # col
    for i in range(n):
        if board[i][col] == 1:
            return False
    # diagonals
    i, j = row-1, col-1
    while i >= 0 and j >= 0:
        if board[i][j] == 1:
            return False
        i -= 1; j -= 1
    i, j = row+1, col+1
    while i < n and j < n:
        if board[i][j] == 1:
            return False
        i += 1; j += 1
    i, j = row-1, col+1
    while i >= 0 and j < n:
        if board[i][j] == 1:
            return False
        i -= 1; j += 1
    i, j = row+1, col-1
    while i < n and j >= 0:
        if board[i][j] == 1:
            return False
        i += 1; j -= 1
    # exactly one queen per color
    color_here = region[row][col]
    for r in range(n):
        for c in range(n):
            if board[r][c] == 1 and region[r][c] == color_here:
                return False
    return True

def solve_with_colors(region):
    """
    Backtracking by rows with constraints:
      - one per row/column/diagonal
      - exactly one per color (color IDs are 0..n-1)
    Returns a solved board or None.
    """
    n = len(region)
    board = [[0]*n for _ in range(n)]
    cols = [False]*n
    diag1 = [False]*(2*n-1)  # r - c + (n-1)
    diag2 = [False]*(2*n-1)  # r + c
    colors_used = [False]*n

    def backtrack(r):
        if r == n:
            return True
        # try columns; prefer columns whose colors unused to reduce dead-ends
        cols_order = list(range(n))
        random.shuffle(cols_order)
        for c in cols_order:
            color = region[r][c]
            d1 = r - c + (n-1)
            d2 = r + c
            if not cols[c] and not diag1[d1] and not diag2[d2] and not colors_used[color]:
                board[r][c] = 1
                cols[c] = diag1[d1] = diag2[d2] = colors_used[color] = True
                if backtrack(r+1):
                    return True
                board[r][c] = 0
                cols[c] = diag1[d1] = diag2[d2] = False
                colors_used[color] = False
        return False

    return board if backtrack(0) else None

def is_partition_solvable(region, max_tries=3):
    """Try a few randomized backtracks to find at least one solution."""
    for _ in range(max_tries):
        sol = solve_with_colors(region)
        if sol is not None:
            return True
    return False

def generate_valid_regions(n, max_attempts=200):
    """Keep generating contiguous blobs until we find a partition with ≥1 solution."""
    for _ in range(max_attempts):
        reg = generate_blob_regions(n)
        if is_partition_solvable(reg, max_tries=3):
            return reg
    return None  # give up (very unlikely for typical N)

# -------------------------------
# GUI
# -------------------------------

class NQueensGUI:
    def __init__(self, n=8):
        self.n = n
        self.palette = generate_palette(n)
        self.region = generate_valid_regions(n)
        if self.region is None:
            raise RuntimeError("Failed to generate a solvable contiguous partition.")
        self.board = [[0]*n for _ in range(n)]  # 0=empty, 1=queen

        self.root = tk.Tk()
        self.root.title(f"{n}-Queens – contiguous blobs (guaranteed ≥1 solution)")

        self.cell_size = max(24, min(72, 720 // max(1, n)))
        self.canvas = tk.Canvas(self.root, width=n*self.cell_size, height=n*self.cell_size)
        self.canvas.pack(pady=6)

        ctrl = tk.Frame(self.root)
        ctrl.pack(pady=4)
        tk.Button(ctrl, text="Check Solution", command=self.check_solution).pack(side=tk.LEFT, padx=4)
        tk.Button(ctrl, text="Auto Solve", command=self.auto_solve).pack(side=tk.LEFT, padx=4)
        tk.Button(ctrl, text="Reset", command=self.reset).pack(side=tk.LEFT, padx=4)
        tk.Button(ctrl, text="New Regions", command=self.new_regions).pack(side=tk.LEFT, padx=4)

        self.canvas.bind("<Button-1>", self.on_click)
        self.draw_board()
        self.root.mainloop()

    # --- Drawing ---
    def draw_board(self):
        self.canvas.delete("all")
        for i in range(self.n):
            for j in range(self.n):
                color_id = self.region[i][j]
                fill = self.palette[color_id % len(self.palette)]
                x1, y1 = j*self.cell_size, i*self.cell_size
                x2, y2 = x1 + self.cell_size, y1 + self.cell_size
                self.canvas.create_rectangle(x1, y1, x2, y2, outline="black", width=1, fill=fill)
                if self.board[i][j] == 1:
                    self.draw_queen(i, j)

    def draw_queen(self, row, col):
        x = col*self.cell_size + self.cell_size//2
        y = row*self.cell_size + self.cell_size//2
        fsize = max(12, int(self.cell_size * 0.55))
        self.canvas.create_text(x, y, text="♛", font=("Arial", fsize), fill="black")

    # --- Interaction ---
    def on_click(self, event):
        col = event.x // self.cell_size
        row = event.y // self.cell_size
        if not (0 <= row < self.n and 0 <= col < self.n):
            return
        if self.board[row][col] == 1:
            self.board[row][col] = 0
            self.draw_board()
            return
        if can_place_manual(self.board, row, col, self.n, self.region):
            self.board[row][col] = 1
            self.draw_board()
        else:
            messagebox.showwarning(
                "Invalid move",
                "Placing a queen here violates a rule:\n"
                "- One per row\n- One per column\n- One per diagonal\n- EXACTLY one per color"
            )

    def check_solution(self):
        n = self.n
        placed = sum(sum(r) for r in self.board)
        if placed != n:
            messagebox.showerror("Not valid", f"You must place exactly {n} queens (currently {placed}).")
            return
        for i in range(n):
            if sum(self.board[i]) != 1:
                messagebox.showerror("Not valid", "A row has 0 or >1 queens.")
                return
        for j in range(n):
            if sum(self.board[i][j] for i in range(n)) != 1:
                messagebox.showerror("Not valid", "A column has 0 or >1 queens.")
                return
        qpos = [(i, self.board[i].index(1)) for i in range(n)]
        d1 = set(); d2 = set(); colors_seen = set()
        for r, c in qpos:
            a = r - c; b = r + c; colid = self.region[r][c]
            if a in d1 or b in d2:
                messagebox.showerror("Not valid", "Queens share a diagonal.")
                return
            if colid in colors_seen:
                messagebox.showerror("Not valid", "A color region has more than one queen.")
                return
            d1.add(a); d2.add(b); colors_seen.add(colid)
        if len(colors_seen) != n:
            messagebox.showerror("Not valid", "Not every color has a queen.")
            return
        messagebox.showinfo("Success", "Valid solution! (Rows, Columns, Diagonals, and Colors)")

    def auto_solve(self):
        sol = solve_with_colors(self.region)
        if sol is None:
            messagebox.showerror("No solution", "Unexpected: these regions had no solution.")
            return
        self.board = sol
        self.draw_board()
        messagebox.showinfo("Solved", "Auto-solved a valid arrangement!")

    def reset(self):
        self.board = [[0]*self.n for _ in range(self.n)]
        self.draw_board()

    def new_regions(self):
        newreg = generate_valid_regions(self.n)
        if newreg is None:
            messagebox.showerror("Error", "Couldn't generate a solvable partition. Try again.")
            return
        self.region = newreg
        self.reset()

# -------------------------------
# Run
# -------------------------------

if __name__ == "__main__":
    NQueensGUI(n=8)
