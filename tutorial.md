# Learning about syndrome gadgets with GadgetLab

This tutorial works through the mathematics behind GadgetLab one object at a time: the code,
its geometric picture, decoding, constrained decoding, cost tables, specs and searches. Each step
comes with a short exercise in the tool. By the end you will be able to set up a gadget search,
read its results, and see in the drawn errors what makes a syndrome pattern behave like a wire
or a copy.

Every number quoted below was computed by GadgetLab, and you can reproduce each one. The worked
examples are saved as session files in `examples/`, so you can open any of them instead of
clicking it together.

What you need: linear algebra over $\mathbb{F}_2$, basic graph theory, and the idea of a CSS
code (X-checks and Z-checks that commute). Nothing else is assumed.

## Contents

1. [Setup](#1-setup)
2. [The code](#2-the-code)
3. [The picture: layers, gaps and faces](#3-the-picture-layers-gaps-and-faces)
4. [Minimum-weight decoding](#4-minimum-weight-decoding)
5. [Constrained decoding: terminals and patches](#5-constrained-decoding-terminals-and-patches)
6. [Cost tables](#6-cost-tables)
7. [Specs: what makes a gadget](#7-specs-what-makes-a-gadget)
8. [Searching for gadgets](#8-searching-for-gadgets)
9. [Local gadgets and what "closed" means](#9-local-gadgets-and-what-closed-means)
10. [A workflow, and the command line](#10-a-workflow-and-the-command-line)
11. [Exercises](#11-exercises)
12. [Reference](#12-reference)

---

## 1. Setup

Install the package as described in [README.md](README.md), then start the engine:

```sh
. .venv/bin/activate
gadgetlab serve
```

Open <http://127.0.0.1:8765>. The viewer has five areas:

- **Top bar:** the code: its name and L, then n, N, checks, k, girth, dX and dZ.
- **Left column:** the graph and number of layers, session buttons, the four modes, and the
  syndrome, region, terminal and patch tools.
- **Canvas:** one panel per layer and per gap (section 3).
- **Right column:** Decode, Cost table, Spec and Search.
- **Legend:** along the bottom, it explains every mark.

Hover over any mark to see its address, for example `Z-check e = 3 (1–2), t = 0`. The four modes
decide what a click does: **1 Syndrome**, **2 Terminals**, **3 Region** and **4 Patch** (the
keys 1 to 4 switch between them). Ctrl+Z undoes the last edit.

To load a session file, press **Open session** and pick it. Results stored in the file are shown
as current when they belong to the session's inputs, and as stale otherwise.

---

## 2. The code

### 2.1 The graph and its cycle code

Let $G$ be a connected cubic graph with vertices $0, \dots, n-1$. Every vertex has degree 3, so
there are $n_E = 3n/2$ edges. GadgetLab stores each edge as $(u, w)$ with $u < w$ and sorts the
list; an edge's index $e$ is its position in that list. For $K_4$ this gives

| e | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| edge | 0–1 | 0–2 | 0–3 | 1–2 | 1–3 | 2–3 |

The incidence matrix $H_1 \in \mathbb{F}_2^{n \times n_E}$ has $H_1[v, e] = 1$ exactly when $v$
is an endpoint of $e$. A vector $c \in \mathbb{F}_2^{n_E}$ is a set of edges, and $H_1 c = 0$
says every vertex meets an even number of them. So $\ker H_1$ is the **cycle space** of $G$:
the edge sets that are edge-disjoint unions of cycles. For a connected graph, $\operatorname{rank} H_1 = n - 1$:
every column has two ones, so the rows sum to zero, and a spanning tree gives $n - 1$
independent columns. Hence

$$
\dim \ker H_1 = n_E - n + 1 = \tfrac{n}{2} + 1, \qquad \dim \ker H_1^{\top} = 1 .
$$

The second number holds because $H_1^\top y = 0$ forces $y_u = y_w$ on every edge, so $y$ is
constant on the connected graph.

### 2.2 The repetition code

$R_L \in \mathbb{F}_2^{(L-1) \times L}$ has ones at $(i, i)$ and $(i, i+1)$. Its kernel is
$\{0^L, 1^L\}$ (dimension 1), and $\ker R_L^\top = \{0\}$.

### 2.3 The hypergraph product

With $\otimes$ the Kronecker product, GadgetLab builds

$$
\begin{aligned}
H_X &= \bigl[\, H_1 \otimes I_L \;\big|\; I_n \otimes R_L^{\top} \,\bigr], \\
H_Z &= \bigl[\, I_{n_E} \otimes R_L \;\big|\; H_1^{\top} \otimes I_{L-1} \,\bigr].
\end{aligned}
$$

The X-checks commute with the Z-checks. By the mixed-product rule
$(A \otimes B)(C \otimes D) = AC \otimes BD$,

$$
\begin{aligned}
H_X H_Z^{\top} &= (H_1 \otimes I_L)(I_{n_E} \otimes R_L^{\top}) + (I_n \otimes R_L^{\top})(H_1 \otimes I_{L-1}) \\
&= H_1 \otimes R_L^{\top} + H_1 \otimes R_L^{\top} = 0 \pmod 2 .
\end{aligned}
$$

The sizes are:

| | count |
|---|---|
| qubits (columns) | $N = n_E L + n (L-1)$ |
| Z-checks (rows of $H_Z$) | $m = n_E (L-1)$ |
| X-checks (rows of $H_X$) | $n L$ |

GadgetLab decodes X errors, so everything from section 4 on uses $H_Z$ only: an error is a 0/1
vector $x$ over the qubits, and its syndrome is $\sigma = H_Z x \bmod 2$.

### 2.4 Code parameters in the top bar

The number of logical qubits is $k = N - \operatorname{rank} H_X - \operatorname{rank} H_Z$, and
GadgetLab computes it from the two ranks. The standard count for hypergraph products,
$k = \dim\ker H_1 \cdot \dim\ker R_L + \dim\ker H_1^\top \cdot \dim\ker R_L^\top$, gives
$(\tfrac{n}{2}+1)\cdot 1 + 1 \cdot 0 = \tfrac{n}{2} + 1$. You can watch this in the tool.

**Try it.** Choose Graph K4 and Layers L 3, then press **New session**. The top bar reads
n 4, N 26, checks 12, k 3, girth 3. Press **Compute distances**: dX 3 and dZ 3. Build Petersen
with L = 3: n 10, N 65, checks 30, k 6, girth 5. Now the distance button is disabled, because
distances are searched exhaustively only up to weight 4, and only when N ≤ 60.

The distances are

$$
\begin{aligned}
d_X &= \min\{\,|x| : H_Z x = 0,\ x \notin \operatorname{row}(H_X)\,\}, \\
d_Z &= \min\{\,|z| : H_X z = 0,\ z \notin \operatorname{row}(H_Z)\,\}.
\end{aligned}
$$

Two operators in the picture of section 3 show where small values come from:

- **One edge in every layer**, $(E, e, 0), \dots, (E, e, L-1)$. It has weight $L$, and every
  face $(e, t)$ contains exactly two of its edges, so $H_Z x = 0$. For every preset and every
  L the engine accepts, it is not a sum of X-checks, so it is an X-type logical operator. For
  K4 the tool reports dX = 2, 3, 4 for L = 2, 3, 4.
- **A cycle of $G$ drawn in a single layer.** It has the cycle's length as weight and meets
  every X-check in an even number of edges, so $H_X z = 0$. A shortest one has length equal to
  the girth. For K4 the tool reports dZ = 3 for every L it computes (2 to 6), and dZ > 4 for
  Petersen (girth 5) at L = 2.

Section 6.2 uses a cost table to show the first of these directly.

---

## 3. The picture: layers, gaps and faces

### 3.1 The product graph

Take $L$ copies of $G$, called layers $0$ to $L-1$. Join copy $t$ to copy $t+1$ by a vertical
edge at every vertex. The result is the product graph $G \,\square\, P_L$, and it carries the
code:

| object | address | index | in the picture |
|---|---|---|---|
| edge-qubit | $(E, e, l)$ | $e L + l$ | edge $e$ in layer $l$ |
| vertex-qubit | $(V, v, t)$ | $n_E L + v(L-1) + t$ | the vertical edge at $v$ between layers $t$ and $t+1$ |
| Z-check | $(e, t)$ | $e(L-1) + t$ | the square face between edge $e$ in layers $t$ and $t+1$ |

**So an error is a set of edges of $G \,\square\, P_L$, and its syndrome is the set of square
faces with an odd number of error edges on their boundary.** The face $(e, t)$ for
$e = (u, w)$ is bounded by exactly four qubits:
$(E, e, t)$, $(E, e, t+1)$, $(V, u, t)$ and $(V, w, t)$.

The canvas draws this one slice at a time. Panel **Layer l** draws $G$ with its edges as the
edge-qubits of layer $l$. Panel **Gap t** draws $G$ again, but now each edge is the face
$(e, t)$ and each vertex is the vertical edge $(V, v, t)$. A fired check is a thick magenta
line; an error is a thick teal edge in a layer panel or a large teal dot in a gap panel.

### 3.2 What a single error fires

Reading off the four-qubit faces gives the column of $H_Z$ for each qubit:

| qubit | fires | column weight |
|---|---|---|
| $(E, e, l)$, $0 < l < L-1$ | $(e, l-1)$ and $(e, l)$ | 2 |
| $(E, e, 0)$ or $(E, e, L-1)$ | only $(e, 0)$, or only $(e, L-2)$ | 1 |
| $(V, v, t)$ | $(e, t)$ for the three edges $e$ at $v$ | 3 |

A helpful way to hold this in your head is a **hypergraph on the faces**. Each face is a node.
An edge-qubit is an ordinary edge joining the face below it to the face above it; in the top
and bottom layers it is a half-edge that ends on the boundary. A vertex-qubit is a hyperedge
joining the three faces around its vertex. Decoding a syndrome then means choosing the fewest
of these edges and hyperedges so that the fired faces are exactly the ones touched an odd
number of times.

**Try it** (K4, L = 3). In Syndrome mode, click the faces $(0,0)$ and $(0,1)$, the edge 0–1 in
Gap 0 and in Gap 1, then press **Decode**. The result is optimal, cost 1, and the witness is
the single edge-qubit $(E, 0, 1)$, drawn in Layer 1. Press **Clear syndrome**, click the three
faces at vertex 0 in Gap 0 (edges 0–1, 0–2, 0–3), and decode again: cost 1, witness
$(V, 0, 0)$, the large teal dot at vertex 0 in Gap 0.

After each decode the viewer recomputes the witness's syndrome itself, from the supports above,
and prints **Syndrome verified**. That check is independent of the solver.

---

## 4. Minimum-weight decoding

### 4.1 The problem

$$
\mathrm{MW}(\sigma) = \min\{\, |x| \;:\; H_Z x = \sigma \,\}.
$$

For this code family $H_Z$ has full row rank $m$ (the spec states it, and the tests check it
for every preset), so
every $\sigma \in \mathbb{F}_2^m$ is a syndrome and $\mathrm{MW}(\sigma)$ is always finite.

GadgetLab solves the problem exactly. The default solver writes each parity condition as an
integer linear constraint with a slack variable,

$$
\sum_{j \in \operatorname{supp}(i)} x_j - 2 y_i = \sigma_i, \qquad x_j \in \{0,1\},\ y_i \in \mathbb{Z}_{\ge 0},
$$

and minimises $\sum_j x_j$ with CBC. A second solver, `brute`, writes the solutions as
$x_0 + \ker$ and enumerates the kernel when its dimension is at most 20. The test suite checks
that the two agree exactly. Several errors can reach the minimum. The tool shows one of them,
and another solver may show a different one of the same weight.

### 4.2 Strings that end on the boundary

Only the top and bottom layers have edge-qubits that fire a single face, so no single qubit
fires a lone face in a middle gap. A string of edge-qubits that runs to one of those layers
does. The string
$(E, e, 0), \dots, (E, e, t)$ fires only $(e, t)$, since every interior face is hit twice; so
does $(E, e, t+1), \dots, (E, e, L-1)$ going up.

**Try it** (K4, L = 4). Fire the single face $(0, 1)$ in Gap 1 and decode: cost 2. The witness
is a vertical string, for example $(E, 0, 2), (E, 0, 3)$ up to the top layer; the string
$(E,0,0), (E,0,1)$ down to the bottom layer costs the same. Fire $(0, 0)$ or $(0, 2)$ instead
and the cost is 1, because those gaps touch a boundary layer.

### 4.3 Strings that meet at a vertex

**Try it** (K4, L = 4). Fire $(0, 1)$ and $(3, 1)$, the faces of edges 0–1 and 1–2 in Gap 1.
Two separate strings would cost 2 + 2. The decoder returns cost 3, for example with witness
$(V, 1, 1), (E, 4, 2), (E, 4, 3)$:

- The vertical edge at vertex 1 fires the three faces around vertex 1, which are
  $(0,1), (3,1), (4,1)$ for the edges 0–1, 1–2 and 1–3.
- The unwanted face $(4, 1)$ is then cancelled by a string from it to the top layer (or, at
  the same cost, to the bottom layer).

---

## 5. Constrained decoding: terminals and patches

Everything else in GadgetLab rests on one primitive. Its inputs are a syndrome $\sigma$, a set
of qubits $Q$ called the patch (all qubits in global scope), and values $\varphi$ for some
terminal qubits $T \subseteq Q$. The problem is

$$
\begin{aligned}
\text{minimise } & \textstyle\sum_{j \in Q} x_j \\
\text{subject to } & x_j = 0 \ (j \notin Q), \quad x_j = \varphi(j) \ (j \in T), \quad (H_Z x)_i = \sigma_i \ (i \in C),
\end{aligned}
$$

where $C$ is the set of enforced checks. In global scope $C$ holds every check. In local_closed
scope it is $C_Q$, the checks whose support meets $Q$. The answer is a cost, or "infeasible"
(shown as "—" and stored as `null`), together with a witness when it is feasible.

### 5.1 Fixing a qubit can raise the cost

**Try it** (K4, L = 3). Fire $(0,0)$ and $(0,1)$ as in section 3. In Terminals mode, click the
edge 0–1 in Layer 1, which makes $(E, 0, 1)$ terminal 1 (a badge "1" appears). Press
**Compute table**:

| a | cost | witness |
|---|---|---|
| 0 | 2 | $(E,0,0), (E,0,2)$ |
| 1 | 1 | $(E,0,1)$ |

With the middle edge forbidden, each fired face needs its own boundary edge.

### 5.2 Parity can make a row infeasible

If every qubit of a fired face is a terminal, the face's parity is decided by the terminal
values alone. **Try it** (K4, L = 3). Fire only $(0, 1)$, then make its four qubits terminals,
in this order: $(E,0,1)$, $(E,0,2)$, $(V,0,1)$, $(V,1,1)$. The table has 16 rows. Every row
with an even number of ones is "—" (infeasible), because then the face sees an even number of
errors. The odd rows are feasible, for example 0100 costs 1, 1000 costs 2, and 0010 and 0001
cost 3.

### 5.3 Patches: local_closed scope

A patch is $Q(U, l_0, l_1)$, built from a vertex set $U$ and a layer range. It holds every
edge-qubit $(E, e, l)$ with $e$ touching $U$ and $l_0 \le l \le l_1$, and every vertex-qubit
$(V, v, t)$ with $v \in U$ and $l_0 \le t \le l_1 - 1$. Local_closed scope means "assume no
error outside the patch". Errors outside $Q$ are forced to 0, and only the checks in $C_Q$ are
enforced.

**Try it** (K4, L = 3). In Patch mode, click vertex 0; set **layers** 1 to 1; choose scope
**local_closed**. The patch is $Q = \{(E,0,1), (E,1,1), (E,2,1)\}$, and $C_Q$ is the six faces
$(e, 0), (e, 1)$ for $e = 0, 1, 2$; hover over a face to see whether it is enforced. Now:

- Fire only $(0,0)$ and decode: **infeasible**. The only patch qubit on that face is
  $(E,0,1)$, and it fires $(0,1)$ as well.
- Add $(0,1)$: cost 1, witness $(E,0,1)$.
- Also fire $(5,0)$, the face of edge 2–3: still cost 1, because that face is not enforced.

---

## 6. Cost tables

### 6.1 Definition and invariants

Choose terminals $T = (t_1, \dots, t_k)$ with $1 \le k \le 6$. For every bit string
$a = a_1 \dots a_k$, sorted with all zeros first, let

$$
\mathrm{Cost}(a) = \text{the constrained cost with } \varphi(t_i) = a_i,
$$

or $\infty$ if that problem is infeasible. GadgetLab also computes the cost with no terminal
fixed and checks two invariants on every table. A violation is an error, never a warning.

1. $\min_a \mathrm{Cost}(a)$ equals the unconstrained cost. Every feasible error gives its
   terminals some values $a$, so the unconstrained feasible set is the disjoint union of the
   row feasible sets.
2. Let $a^*$ be the terminal values of the unconstrained witness. Then
   $\mathrm{Cost}(a^*)$ equals the unconstrained cost, because that witness is feasible for
   row $a^*$.

**Try it.** Open `examples/petersen_L3.session.json`. Its decode cost is 3, and its cost table
reads 00: 3, 01: 5, 10: 4, 11: 4. The minimum, 3, is the decode cost. The decode witness
$(E,0,0), (E,3,2), (E,11,0)$ sets both terminals, $(E,3,1)$ and $(V,1,0)$, to 0, and row 00
costs 3.

### 6.2 A table with no syndrome: stabilizers and logicals

With $\sigma = 0$, row $a$ asks for the lightest zero-syndrome error with the given terminal
values. **Try it** (K4, L = 3). Press **Clear syndrome**, then make the column of edge 0–1
terminals: $(E,0,0)$, $(E,0,1)$, $(E,0,2)$. **Compute table**. Where several errors tie, your
copy of the tool may draw a different witness of the same weight; the costs do not change.

| a | cost | the witness the tool draws |
|---|---|---|
| 000 | 0 | no error |
| 100 | 4 | the X-check at vertex 0 in layer 0: $(E,0,0),(E,1,0),(E,2,0),(V,0,0)$ |
| 010 | 5 | the X-check at vertex 1 in layer 1 (a middle layer, so weight 5) |
| 001 | 4 | the X-check at vertex 1 in layer 2 |
| 111 | 3 | the column itself: $(E,0,0),(E,0,1),(E,0,2)$ |
| 011 | 5 | the column plus the X-check at vertex 0 in layer 0; their common qubit cancels |
| 101 | 6 | the column plus the X-check at vertex 1 in layer 1 |
| 110 | 5 | the column plus the X-check at vertex 1 in layer 2 |

Zero-syndrome errors come in two kinds, and the table shows both. Products of X-checks
(stabilizers) are trivial. The column of weight $L = 3$ is a logical operator, and the tool's
dX = 3 for this code says nothing lighter is.

---

## 7. Specs: what makes a gadget

### 7.1 The idea

The research project behind GadgetLab looks for **gadgets**: small syndrome patterns whose
decoding cost behaves like a logic component, for example a wire that carries a bit, or a
branch that copies it. "Behaves like" is made precise by asking for a cost table of a given
shape.

### 7.2 The spec language

A spec gives each row $a$ exactly one entry, where $k_a$ is a non-negative integer:

| entry | meaning |
|---|---|
| `["eq", k_a]` | $\mathrm{Cost}(a) = c + k_a$ |
| `["ge", k_a]` | $\mathrm{Cost}(a) \ge c + k_a$, and $\infty$ satisfies it |
| `["any"]` | no condition |

At least one row must be eq. A table **matches** the spec when a single integer $c$ makes every
entry true. That $c$ is read off the eq rows as $\mathrm{Cost}(a) - k_a$; every eq row must be
finite and give the same value. The **margin** measures how far the ge rows clear their bound:

$$
\text{margin} = \min_{a \,:\, \text{ge},\ \mathrm{Cost}(a) < \infty} \bigl(\mathrm{Cost}(a) - c - k_a\bigr),
$$

and it is $\infty$ (stored as `null`) when no ge row is finite. When the eq rows disagree, the
table does not match and $c$ is undefined. When they agree but a ge row falls short, the table
does not match, and the margin shows the shortfall as a negative number.

The presets are:

- **wire** (2 terminals): 00 eq 0, 11 eq 0, 01 ge 1, 10 ge 1. The cheapest explanations of the
  syndrome give both terminals the same value, either way at the same cost $c$. Any explanation
  where they differ costs at least $c + 1$, and the margin is how much more than that.
- **copy3** (3 terminals): 000 eq 0, 111 eq 0, every other row ge 1. All three terminals agree
  in the cheapest explanations.
- **custom**: every row "any" until you edit it. Changing the number of terminals resets the
  spec to custom, and the viewer says so.

### 7.3 A wire on K4

Open `examples/k4_L2.session.json` (K4, L = 2). Its terminals are the vertical edges at
vertices 0 and 2, $(V,0,0)$ and $(V,2,0)$, and its spec is wire. The stored search found two
patterns. Click result **1** in the Search section: the canvas shows "From search result 1" and
fires the four faces around the cycle 0–1–2–3, and the cost table loads:

| a | cost | witness |
|---|---|---|
| 00 | 2 | $(V,1,0), (V,3,0)$ |
| 01 | 4 | |
| 10 | 4 | |
| 11 | 2 | $(V,0,0), (V,2,0)$ |

The table reads **Matches the spec: c = 2, margin 1**. Press **Show** on rows 00 and 11 to see
the two cheapest explanations:

- Row 11 draws the vertical edges at 0 and 2. Each fires the three faces around its vertex.
  The face of edge 0–2 is fired twice and cancels, leaving exactly the four faces of the cycle.
- Row 00 does the same with the other two corners, 1 and 3, using the face of edge 1–3.
- The terminals sit on opposite corners of the cycle, so they are both on in one explanation
  and both off in the other. Mixed values cost $4 = c + 1 + 1$: the bound plus margin 1.

### 7.4 Same syndrome, other terminals

Open `examples/k4_L2_opposite.session.json`. It has the same four fired faces, but the
terminals are at adjacent corners, 0 and 1, and the spec is custom: 01 and 10 eq 0, 00 and 11
ge 1. The table is 00: 4, 01: 2, 10: 2, 11: 6, and it **matches with c = 2, margin 1**. The
cheap explanations are the same two pairs of corners as before. Each pair contains exactly one
of the terminals, so in each cheapest explanation the terminals disagree.

### 7.5 A copy on the prism

Open `examples/prism_L2_copy3.session.json`. The prism has an outer triangle 0, 1, 2, an inner
triangle 3, 4, 5, and three spokes 0–3, 1–4, 2–5. The fired faces are the three spokes, the
terminals are the vertical edges at 0, 1 and 2, and the spec is copy3. The table:

| a | cost | witness |
|---|---|---|
| 000 | 3 | for example $(V,3,0), (V,4,0), (V,5,0)$, the inner triangle |
| 111 | 3 | $(V,0,0), (V,1,0), (V,2,0)$, the outer triangle |
| any other row | 5 | |

So the table matches copy3 with c = 3 and margin 1. Around a triangle every triangle edge is
shared by two chosen vertices and cancels, which leaves only the spokes. Row 000 has other
cheapest explanations too: one edge-qubit on each spoke also costs 3.

---

## 8. Searching for gadgets

### 8.1 What a search does

A search fixes a base syndrome $\sigma_0$, terminals, a spec and a scope. It also takes a
**region** $R$: an ordered list of checks whose values it varies. A candidate is a pattern
$p \in \{0,1\}^{|R|}$, and its syndrome is $\sigma_0$ with the entries on $R$ replaced by $p$.

- If $2^{|R|} \le B$ for the budget $B$ (default 4096), every pattern is tried; the run is
  **exhaustive**. Otherwise $B$ distinct patterns are drawn uniformly with the seed.
- Patterns with more than **max weight** ones are skipped and counted as filtered.
- Evaluation is lazy. The eq rows are decoded first and fix $c$, rejecting at the first
  disagreement. Then the ge rows are decoded, rejecting at the first violation. Only matches
  get their remaining rows filled in.
- Matches are ranked by margin (largest first, $\infty$ first of all), then $c$ (smallest
  first), then $|p|$ (fewest ones first), then $p$. The top $k$ are kept with full tables and a
  witness per row.

### 8.2 Running one in the viewer

1. Set terminals (Terminals mode) and a spec (the Spec section).
2. In Region mode (key 3), click the faces to vary. Each gets a dashed halo, and the Search
   section reports the region size and whether the run will be exhaustive or sampled, for
   example "Region: 10 checks. 2^10 = 1,024 patterns: exhaustive."
3. Set budget, seed, top k, max weight and workers, or leave the defaults, and press
   **Run search**. The progress bar and the counts update live; **Cancel search** stops it.
4. When it finishes, click a row of the results list (rank, c, margin, |p|, pattern) to load
   its syndrome and cost table. Then use **Show** on the table rows to see why it matches.

The summary counts add up: candidates = filtered + infeasible syndromes + rejected + matched.
In global scope the infeasible-syndrome count is 0 for this family, because $H_Z$ has full row
rank.

### 8.3 A larger example

Open `examples/petersen_L3.session.json`. It is Petersen with L = 3, three fired faces, the
terminals $(E,3,1)$ and $(V,1,0)$, the wire spec, and a region of 10 checks. The stored search
tried all 1024 patterns and found 92 matches. The list is sorted by margin and rank 1 has
margin 0, so every match has margin 0. Result 1 adds the two faces of edge 1–6, in Gap 0 and
Gap 1. Its table reads 00: 4, 01: 5, 10: 5, 11: 4, so c = 4 and margin 0. Compare the
witnesses of rows 00 and 11 to see which part of the error each terminal value switches.

Two habits help when reading results. First, a margin-0 wire still meets the spec, but its
mixed rows sit exactly on the bound $c + 1$; a larger margin means the mixed rows cost more.
Second, the shown witness is one optimal error among possibly several, so compare costs, not
pictures, when you judge a pattern.

---

## 9. Local gadgets and what "closed" means

In local_closed scope the search only enforces the checks in $C_Q$. Region checks outside $C_Q$
are dropped with a warning. They keep their base value in the stored pattern but are not
varied.

Open `examples/prism_L2_local.session.json`. The patch is $U = \{0,1,2,3\}$, which is
ball(0, 1) (the **Fill patch with ball** helper with vertex 0 and radius 1), with layers 0 to 1.
That gives $Q$ = 20 qubits and $C_Q$ = 8 faces: every Gap 0 face except the one of edge 4–5.
The terminals are the vertical edges at 0 and 1, and the spec is wire. The stored search tried
256 patterns and found 8 matches, all with margin 1.

The session's own syndrome is result 3: the faces of edges 1–4, 2–5, 3–4 and 3–5. Its local
table matches the wire with c = 4:

| a | local cost | witness |
|---|---|---|
| 00 | 4 | for example $(E,4,0), (E,5,1), (E,6,1), (E,7,1)$, one edge-qubit per fired face |
| 01, 10 | 6 | |
| 11 | 4 | $(V,0,0), (V,1,0), (V,2,0), (V,3,0)$ |

Now switch scope to **global**. The table and the decode turn stale, because scope is one of
their inputs. Press **Decode**: cost 2, with witness $(V,4,0), (V,5,0)$. Those are the
vertical edges at vertices 4 and 5, both outside the patch. **Compute table** gives
00: 2, 01: 6, 10: 6, 11: 4: the eq rows disagree, so in global scope this pattern is not a
wire.

That is what "closed" means. A local_closed result holds under the assumption that no error
appears outside the patch, and it can fail once errors there are allowed. Of the 8 local
matches in this example, 7 still match in global scope; result 3 is the exception. GadgetLab v1
does not test robustness against errors just outside a patch; that is planned for v2.

---

## 10. A workflow, and the command line

A loop that works well:

1. Start small: K4, K3,3 or the prism with L = 2 or 3, and a region of up to about 10 checks,
   so searches are exhaustive and take seconds.
2. Pick terminals and a spec, search, and load the best results.
3. For each result, press **Show** on every row and compare the witnesses. The difference
   between the eq rows and the ge rows is the mechanism.
4. Vary one thing at a time: where the terminals sit, L, the scope and patch, the spec's
   offsets. Undo (Ctrl+Z) brings back earlier inputs, and a stored result becomes current again
   when its inputs return.
5. Save the session (**Save session**) and reproduce it from the command line.

The command line computes everything the viewer can, and writes results into the session file.
Start with a fresh session:

```sh
gadgetlab build --preset prism -L 2 --out prism.json
```

To script the inputs, edit the session JSON. Addresses are written as in §5 of the spec, for
example `{"e": 3, "t": 0}` for a check and `{"q": "V", "v": 2, "t": 0}` for a vertex-qubit.
This sets up the copy3 example of section 7.5:

```python
import json
s = json.load(open("prism.json"))
s["syndrome"] = [{"e": e, "t": 0} for e in (2, 4, 5)]           # the three spokes
s["terminals"] = [{"q": "V", "v": v, "t": 0} for v in (0, 1, 2)]
s["spec"] = {"preset": "copy3", "rows": {f"{i:03b}": (["eq", 0] if i in (0, 7) else ["ge", 1]) for i in range(8)}}
s["region"] = [{"e": e, "t": 0} for e in range(9)]
json.dump(s, open("prism.json", "w"), indent=2)
```

Then compute:

```sh
gadgetlab decode prism.json                  # cost 3
gadgetlab table prism.json                   # 000 and 111 cost 3, the rest 5: matches, c = 3, margin 1
gadgetlab search prism.json --workers 3      # 512 candidates, 4 matched
gadgetlab search prism.json --check          # recompute and compare with the stored result
```

Searches are deterministic for a fixed seed, and the CLI writes the same result JSON as the
viewer's engine, apart from timing fields.

---

## 11. Exercises

Answers are below each one; try first, then check with the tool.

**1.** On K4 with L = 3, which faces does the vertical edge $(V, 2, 1)$ fire? Fire them and
decode.

<details><summary>Answer</summary>

The faces of the three edges at vertex 2 (0–2, 1–2, 2–3) in Gap 1: $(1,1), (3,1), (5,1)$.
Decoding them returns cost 1 with witness $(V,2,1)$.
</details>

**2.** On K4 with L = 5, fire a single face $(0, t)$ for $t = 0, 1, 2, 3$ in turn. Predict each
cost from section 4.2, then decode.

<details><summary>Answer</summary>

Costs 1, 2, 2, 1. The witnesses are $(E,0,0)$; then $(E,0,0),(E,0,1)$; then
$(E,0,3),(E,0,4)$; then $(E,0,4)$: each is a string to the nearer boundary layer.
</details>

**3.** On K4 with L = 4, fire $(0, 0)$ and $(0, 2)$. Is one string between them cheaper than
two strings to the boundary?

<details><summary>Answer</summary>

No, it is not cheaper. The cost is 2, for example with witness $(E,0,0), (E,0,3)$: one boundary
edge for each face. A string between them, $(E,0,1),(E,0,2)$, also costs 2.
</details>

**4.** In `examples/k4_L2_opposite.session.json`, move terminal 2 from vertex 1 to vertex 2:
remove it with × in the Cost table section, then in Terminals mode (key 2) click vertex 2 in
Gap 0. Which spec matches now?

<details><summary>Answer</summary>

Moving a terminal changes their number twice, so the spec resets to custom; choose **wire**.
The table is 00: 2, 01: 4, 10: 4, 11: 2, which matches the wire with c = 2 and margin 1.
This is the situation of section 7.3.
</details>

**5.** In `examples/prism_L2_local.session.json`, take each of the 8 search results in turn:
with scope local_closed, click the result; then switch scope to global and press
**Compute table**. Which ones stop matching? (Switch back to local_closed before loading the
next result. In global scope the search itself is stale, so its results cannot be loaded.)

<details><summary>Answer</summary>

Only result 3, the session's own syndrome (section 9). The other 7 still match the wire in
global scope, with the same c and margin.
</details>

**6.** Check the two invariants of section 6.1 on any table you have computed. Then explain in
one sentence why invariant 1 cannot fail for a correct solver.

<details><summary>Answer</summary>

Every feasible error, constrained or not, gives the terminals some values, so the unconstrained
optimum is the best of the row optima.
</details>

---

## 12. Reference

### Index map

As in §3.5 of the spec:

| object | address | index |
|---|---|---|
| edge-qubit | $(E, e, l)$, $0 \le l \le L-1$ | $eL + l$ |
| vertex-qubit | $(V, v, t)$, $0 \le t \le L-2$ | $n_E L + v(L-1) + t$ |
| Z-check | $(e, t)$, $0 \le t \le L-2$ | $e(L-1) + t$ |
| X-check | $(v, l)$, $0 \le l \le L-1$ | $vL + l$ |

### Viewer marks

| mark | meaning |
|---|---|
| thick magenta edge (gap panel) | check with σ = 1 |
| thin grey edge | check with σ = 0, or an idle edge-qubit (layer panel) |
| thick teal edge (layer panel) | error on an edge-qubit |
| large teal dot (gap panel) | error on a vertex-qubit |
| numbered badge | terminal |
| dashed halo | region check |
| light blue shading, dashed ring | patch qubits, patch vertex |

### Words used in this tutorial

| word | meaning |
|---|---|
| witness | an optimal error found by the solver |
| cost | the weight of the witness inside the patch, or "—" if infeasible |
| row | one assignment $a$ of terminal values |
| c | the common baseline cost of the eq rows |
| margin | how far the ge rows clear their bound; ∞ if none is finite |
| exhaustive | the search tried every pattern on the region |

For the precise definitions, see §3 of [gadgetlab-v1-spec.md](gadgetlab-v1-spec.md); for the
commands and the HTTP API, see [README.md](README.md).
