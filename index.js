// Define physics constants
const gravity = 0.06666666666666665;
const gravity_held = 0.01111111111111111;
const ground_accel = 0.06666666666666665;
const air_accel = 0.04444444444444444;
const drag = 0.9933221725495059; // 0.99^(2/3)
const friction_ground = 0.9459290248857720; // 0.92^(2/3)
const friction_ground_slow = 0.8617738760127536; // 0.80^(2/3)
const friction_wall = 0.9113380468927672; // 0.87^(2/3)
const max_xspeed = 3.333333333333333;
const max_jump_duration = 45;

// Define map variables
let mdata, tile_data, 
    tile_dic, segment_dic, entity_dic, entity_list,
    hor_grid_edge_dic, ver_grid_edge_dic, hor_segment_dic, ver_segment_dic,
    tile_grid_edge_map, tile_segment_map;

// Define ninja
let p;

// Define framerate variables and track time of last frame
let fps = 70;
let interval = 1000 / fps;
let unlimited = false;
let last_frame = Date.now();

// Define canvas
const canvas = document.querySelector("canvas");
const ctx = canvas.getContext("2d");

// Track inputs
let input = { left: 0, right: 0, jump: 0 };
document.addEventListener('keydown', key => {
    if (key.repeat) return;
    switch (key.code) {
        case 'ArrowLeft':  input.left  = 1; break;
        case 'ArrowRight': input.right = 1; break;
        case 'KeyZ':       input.jump  = 1; break;
    }
});
document.addEventListener('keyup', key => {
    switch (key.code) {
        case 'ArrowLeft':  input.left  = 0; break;
        case 'ArrowRight': input.right = 0; break;
        case 'KeyZ':       input.jump  = 0; break;
    }
});

// This class is responsible for updating and storing the positions and velocities of each ninja
// The properties xposlog and yposlog contain all the coordinates used to generate the traces of the replays
class Ninja {
    // Initiate ninja position at spawn point, and initiate other values to their initial state
    constructor(xspawn, yspawn) {
        this.xpos = xspawn;
        this.ypos = yspawn;
        this.xspeed = 0;
        this.yspeed = 0;
        this.applied_gravity = gravity;
        this.applied_friction = friction_ground;
        this.state = 0; // 0: Immobile, 1: Running, 2: Ground sliding, 3: Jumping, 4: Falling, 5: Wall sliding
        this.radius = 10;
        this.hor_input = 0;
        this.hor_input = 0;
        this.jump_input = 0;
        this.jump_input_old = 0;
        this.airborn = true;
        this.walled = false;
        this.jump_duration = 0;
        this.jump_buffer = -1;
        this.floor_buffer = -1;
        this.wall_buffer = -1;
        this.launch_pad_buffer = -1;
        this.poslog = [0, xspawn, yspawn]; // Used for debug
        this.speedlog = [0, 0, 0];
        this.xposlog = [xspawn]; // Used to produce trace
        this.yposlog = [yspawn];
    }

    // Update position and speed by applying drag and gravity before collision phase
    integrate() {
        this.xspeed *= drag;
        this.yspeed *= drag;
        this.yspeed += this.applied_gravity;
        this.xpos_old = this.xpos;
        this.ypos_old = this.ypos;
        this.xpos += this.xspeed;
        this.ypos += this.yspeed;
    }

    // Reset some values used for collision phase
    pre_collision() {
        this.xspeed_old = this.xspeed;
        this.yspeed_old = this.yspeed;
        this.floor_count = 0;
        this.wall_count = 0;
        this.floor_normal_x = 0;
        this.floor_normal_y = 0;
    }

    // Gather all entities in neighbourhood and apply physical collisions if possible
    collide_vs_objects() {
        // TODO
    }

    // Gather all tile segments in neighbourhood and handle collisions with those
    collide_vs_tiles() {
        let dx = this.xpos - this.xpos_old;
        let dy = this.ypos - this.ypos_old;

        // Interpolation routine mainly to prevent from going through walls
        const time = sweep_circle_vs_tiles(this.xpos_old, this.ypos_old, dx, dy, this.radius * 0.5);
        this.xpos = this.xpos_old + time * dx;
        this.ypos = this.ypos_old + time * dy;

        // Find the closest point from the ninja, apply depenetration and update speed. Loop 32 times.
        for (let point = 0; point < 32; point++) {
            const closest_point = get_single_closest_point(this.xpos, this.ypos, this.radius);
            if (!closest_point) break;
            const a = closest_point[0];
            const b = closest_point[1];
            let dx = this.xpos - a;
            let dy = this.ypos - b;
            if (Math.abs(dx) <= 0.0000001) // This is to prevent corner cases. Not in the original code.
                dx = 0;
            if (this.xpos == 50.51197510492316 || this.xpos == 49.23232124849253) // This is to artificially create corner cases at certain common exact spots. lol
                dx = Math.pow(-2, -47);
            if (this.xpos == 49.153536108584795)
                dx = Math.pow(2, -47);
            const dist = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
            if (dist == 0) break;
            const xpos_new = a + this.radius * dx / dist;
            const ypos_new = b + this.radius * dy / dist;
            this.xpos = xpos_new;
            this.ypos = ypos_new;
            const dot_product = this.xspeed * dx + this.yspeed * dy;
            if (dot_product < 0) { // Check if you're moving towards the corner
                const xspeed_new = (this.xspeed * dy - this.yspeed * dx) / Math.pow(dist, 2) * dy;
                const yspeed_new = (this.xspeed * dy - this.yspeed * dx) / Math.pow(dist, 2) * -dx;
                this.xspeed = xspeed_new;
                this.yspeed = yspeed_new;
            }
            if (dy < -0.0001) {
                this.floor_count += 1;
                this.floor_normal_x += dx / dist;
                this.floor_normal_y += dy / dist;
            }
        }
    }

    // Perform logical collisions with entities, check for walled state and calculate floor normals
    post_collision() {
        // Perform LOGICAL collisions between the ninja and nearby entities.
        // Also check if the ninja can interact with the walls of entities when applicable.
        let wall_normal = 0;
        // TODO

        // Check if the ninja can interact with walls from nearby tile segments
        const rad = this.radius + 0.1;
        const neighbour_cells = gather_cells_from_region(this.xpos - rad, this.ypos - rad, this.xpos + rad, this.ypos + rad);
        for (const cell of neighbour_cells) {
            for (const segment of segment_dic[cell[0]][cell[1]]) {
                if (segment.active && segment.type == "linear") {
                    const collision_result = segment.wall_intersecting(this);
                    if (collision_result) wall_normal += collision_result;
                }
            }
        }

        // Calculate combined wall normal
        this.airborn = true;
        this.walled = false;
        if (wall_normal) {
            this.walled = true;
            this.wall_normal = wall_normal / Math.abs(wall_normal);
        }

        // Calculate the combined floor normalized normal vector if the ninja has touched any floor
        if (this.floor_count > 0) {
            this.airborn = false;
            const floor_scalar = Math.sqrt(Math.pow(this.floor_normal_x, 2) + Math.pow(this.floor_normal_y, 2));
            if (floor_scalar == 0) {
                this.floor_normalized_x = 0;
                this.floor_normalized_y = -1;
            } else {
                this.floor_normalized_x = this.floor_normal_x / floor_scalar;
                this.floor_normalized_y = this.floor_normal_y / floor_scalar;
            }
        }
    }

    // Perform floor jump depending on slope angle and direction
    floor_jump() {
        this.jump_buffer = -1;
        this.floor_buffer = -1;
        this.launch_pad_buffer = -1;
        this.state = 3;
        this.applied_gravity = gravity_held;
        let jx, jy;
        if (this.floor_normalized_x == 0) { // Jump from flat ground
            jx = 0;
            jy = -2;
        } else { // Slope jump
            const dx = this.floor_normalized_x;
            const dy = this.floor_normalized_y;
            if (this.xspeed * dx >= 0) { // Moving downhill
                if (this.xspeed * this.hor_input >= 0) {
                    jx = 2 / 3 * dx;
                    jy = 2 * dy;
                } else {
                    jx = 0;
                    jy = -1.4;
                }
            } else { // Moving uphill
                if (this.xspeed * this.hor_input > 0) { // Forward jump
                    jx = 0;
                    jy = -1.4;
                } else {
                    this.xspeed = 0; // Perp jump
                    jx = 2 / 3 * dx;
                    jy = 2 * dy;
                }
            }
        }
        if (this.yspeed > 0) this.yspeed = 0;
        this.xspeed += jx;
        this.yspeed += jy;
        this.xpos += jx;
        this.ypos += jy;
        this.jump_duration = 0;
    }

    // Perform wall jump depending on wall normal and if sliding or not
    wall_jump() {
        let jx, jy;
        if (this.hor_input * this.wall_normal < 0 && this.state == 5) { // Slide wall jump
            jx = 2 / 3;
            jy = -1;
        } else { // Regular wall jump
            jx = 1;
            jy = -1.4;
        }
        this.state = 3;
        this.applied_gravity = gravity_held;
        if (this.xspeed * this.wall_normal < 0) this.xspeed = 0;
        if (this.yspeed > 0) this.yspeed = 0;
        this.xspeed += jx * this.wall_normal;
        this.yspeed += jy;
        this.xpos += jx * this.wall_normal;
        this.ypos += jy;
        this.jump_buffer = -1;
        this.wall_buffer = -1;
        this.launch_pad_buffer = -1;
        this.jump_duration = 0;
    }

    // Perform launch pad jump
    lp_jump() {
        this.floor_buffer = -1;
        this.wall_buffer = -1;
        this.jump_buffer = -1;
        this.launch_pad_buffer = -1;
        let boost_scalar = 2 * Math.abs(self.xlp_boost_normalized) + 2;
        if (boost_scalar == 2) boost_scalar = 1.7; // This was really needed. Thanks Metanet
        this.xspeed += this.xlp_boost_normalized * boost_scalar * 2 / 3;
        this.yspeed += this.ylp_boost_normalized * boost_scalar * 2 / 3;
    }

    // This function handles all of the ninja's actions depending on the inputs and the environment
    think() {
        // Logic to determine if you're starting a new jump
        let new_jump_check;
        if (!this.jump_input)
            new_jump_check = false;
        else
            new_jump_check = this.jump_input_old == 0;
        this.jump_input_old = this.jump_input;

        // Determine if within buffer ranges. If so, increment buffers.
        if (-1 < this.launch_pad_buffer && this.launch_pad_buffer < 3)
            this.launch_pad_buffer += 1;
        else
            this.launch_pad_buffer = -1;
        const in_lp_buffer = -1 < this.launch_pad_buffer && this.launch_pad_buffer < 4;
        if (-1 < this.jump_buffer && this.jump_buffer < 5)
            this.jump_buffer += 1;
        else
            this.jump_buffer = -1;
        const in_jump_buffer = -1 < this.jump_buffer && this.jump_buffer < 5;
        if (-1 < this.wall_buffer && this.wall_buffer < 5)
            this.wall_buffer += 1;
        else
            this.wall_buffer = -1;
        const in_wall_buffer = -1 < this.wall_buffer && this.wall_buffer < 5;
        if (-1 < this.floor_buffer && this.floor_buffer < 5)
            this.floor_buffer += 1;
        else
            this.floor_buffer = -1;
        const in_floor_buffer = -1 < this.floor_buffer && this.floor_buffer < 5;

        // Initiate jump buffer if beginning a new jump and airborn
        if (new_jump_check && this.airborn) this.jump_buffer = 0;
        // Initiate wall buffer if touched a wall this frame
        if (this.walled) this.wall_buffer = 0;
        // Initiate floor buffer if touched a floor this frame
        if (!this.airborn) this.floor_buffer = 0;

        // This block deals with the case where the ninja is touching a floor
        if (!this.airborn) {
            const xspeed_new = this.xspeed + ground_accel * this.hor_input;
            if (Math.abs(xspeed_new) < max_xspeed) this.xspeed = xspeed_new;
            if (this.state > 2) {
                if (this.xspeed * this.hor_input <= 0) {
                    if (this.state == 3) this.applied_gravity = gravity;
                    this.state = 2;
                } else {
                    if (this.state == 3) this.applied_gravity = gravity;
                    this.state = 1;
                }
            }
            if (!in_jump_buffer && !new_jump_check) { // If not jumping
                if (this.state == 2) {
                    const projection = Math.abs(this.yspeed * this.floor_normalized_x - this.xspeed * this.floor_normalized_y);
                    if (this.hor_input * projection * this.xspeed > 0) {
                        this.state = 1;
                        return;
                    }
                    if (projection < 0.1 && this.floor_normalized_x == 0) {
                        this.state = 0;
                        return;
                    }
                    if (this.yspeed < 0 && this.floor_normalized_x != 0) {
                        // Up slope friction formula, very dumb but that's how it is
                        const speed_scalar = Math.sqrt(Math.pow(this.xspeed, 2) + Math.pow(this.yspeed, 2));
                        const fric_force = Math.abs(this.xspeed * (1 - friction_ground) * this.floor_normalized_y);
                        const fric_force2 = speed_scalar - fric_force * Math.pow(this.floor_normalized_y, 2);
                        this.xspeed = this.xspeed / speed_scalar * fric_force2;
                        this.yspeed = this.yspeed / speed_scalar * fric_force2;
                        return;
                    }
                    this.xspeed *= friction_ground;
                    return;
                }
                if (this.state == 1) {
                    const projection = Math.abs(this.yspeed * this.floor_normalized_x - this.xspeed * this.floor_normalized_y);
                    if (this.hor_input * projection * this.xspeed > 0) {
                        if (this.hor_input * this.floor_normalized_x >= 0) return; // If holding inputs in downhill direction or flat ground
                        if (Math.abs(xspeed_new) < max_xspeed) {
                            const boost = ground_accel / 2 * this.hor_input;
                            const xboost = boost * this.floor_normalized_y * this.floor_normalized_y;
                            const yboost = boost * this.floor_normalized_y * -this.floor_normalized_x;
                            this.xspeed += xboost;
                            this.yspeed += yboost;
                        }
                        return;
                    }
                    this.state = 2;
                } else { // If you were in state 0 I guess
                    if (this.hor_input) {
                        this.state = 1;
                        return;
                    }
                    const projection = Math.abs(this.yspeed * this.floor_normalized_x - this.xspeed * this.floor_normalized_y);
                    if (projection < 0.1) {
                        this.xspeed *= friction_ground_slow;
                        return;
                    }
                    this.state = 2;
                }
                return;
            }
            // If you're jumping
            this.floor_jump();
            return;
        } else { // This block deals with the case where the ninja didn't touch a floor
            const xspeed_new = this.xspeed + air_accel * this.hor_input;
            if (Math.abs(xspeed_new) < max_xspeed) this.xspeed = xspeed_new;
            if (this.state < 3) {
                this.state = 4;
                return;
            }
            if (this.state == 3) {
                this.jump_duration += 1;
                if (!this.jump_input || this.jump_duration > max_jump_duration) {
                    this.applied_gravity = gravity;
                    this.state = 4;
                    return;
                }
            }
            if (in_jump_buffer || new_jump_check) {
                if (this.walled || in_wall_buffer) {
                    this.wall_jump();
                    return;
                }
                if (in_floor_buffer) {
                    this.floor_jump();
                    return;
                }
                if (in_lp_buffer && new_jump_check) {
                    this.lp_jump();
                    return;
                }
            }
            if (!this.walled) {
                if (this.state == 5) this.state = 4;
            } else {
                if (this.state == 5) {
                    if (this.hor_input * this.wall_normal <= 0)
                        this.yspeed *= friction_wall;
                    else
                        this.state = 4;
                } else {
                    if (this.yspeed > 0 && this.hor_input * this.wall_normal < 0) {
                        if (this.state == 3) this.applied_gravity = gravity;
                        this.state = 5;
                    }
                }
            }
        }
    }
}

// Contains all the linear segments of tiles and doors that the ninja can interact with
class GridSegmentLinear {
    // Initiate an instance of a linear segment of a tile.
    // Each segment is defined by the coordinates of its two end points.
    constructor(p1, p2) {
        [this.x1, this.y1] = p1;
        [this.x2, this.y2] = p2;
        this.active = true;
        this.type = "linear";
    }

    // Check if the ninja is interesecting with the segment.
    // If so, return the penetration length and the closest point on the segment from the center of the ninja.
    collision_check(xpos, ypos, radius) {
        const px = this.x2 - this.x1;
        const py = this.y2 - this.y1;
        const seg_lensq = Math.pow(px, 2) + Math.pow(py, 2);
        const u = Math.max(0, Math.min(1, ((xpos - this.x1) * px + (ypos - this.y1) * py) / seg_lensq));
        const x = this.x1 + u * px;
        const y = this.y1 + u * py;
        const dist = Math.sqrt(Math.pow(xpos - x, 2) + Math.pow(ypos - y, 2));
        const penetration = dist < 9.9999999 ? (radius - dist) : 0;
        return [[x, y], penetration];
    }

    // Return the time of intersection (as a fraction of a frame) for the closest point in the ninja's path.
    // Return 0 if the ninja is already intersection or 1 if won't intersect within the frame.
    intersect_with_ray(xpos, ypos, dx, dy, radius) {
        const time1 = get_time_of_intersection_circle_vs_circle(xpos, ypos, dx, dy, this.x1, this.y1, radius);
        const time2 = get_time_of_intersection_circle_vs_circle(xpos, ypos, dx, dy, this.x2, this.y2, radius);
        const time3 = get_time_of_intersection_circle_vs_lineseg(xpos, ypos, dx, dy, this.x1, this.y1, this.x2, this.y2, radius);
        return Math.min(time1, time2, time3);
    }

    // If the ninja with an icreased radius of 10.1 is intersecting a wall, return the wall normal
    wall_intersecting(ninja) {
        if (this.x1 == this.x2 && this.y1 <= ninja.ypos && ninja.ypos <= this.y2) {
            if (-(ninja.radius + 0.1) < ninja.xpos - this.x1 && ninja.xpos - this.x1 < 0)
                return -1;
            if (0 < ninja.xpos - this.x1 && ninja.xpos - this.x1 < (ninja.radius + 0.1))
                return 1;
        }
    }
}

// Contains all the circular segments of tiles that the ninja can interract with
class GridSegmentCircular {
    // Initiate an instance of a circular segment of a tile. Each segment is defined by:
    // - The coordinates of its center
    // - A vector indicating which quadrant contains the qurater-circle
    // - A boolean indicating if the tile is convex or concave
    // - The radius of the quarter-circle
    constructor(center, quadrant, convex = true, radius = 24) {
        this.xpos = center[0];
        this.ypos = center[1];
        this.hor = quadrant[0];
        this.ver = quadrant[1];
        this.active = true;
        this.type = "circular";
        this.radius = radius;
        this.convex = convex;
    }

    // Check if the ninja is interesecting with the segment.
    // If so, calculate the penetration length and the closest point on the segment from the center of the ninja.
    collision_check(xpos, ypos, radius) {
        const dx = xpos - this.xpos;
        const dy = ypos - this.ypos;
        let x, y, penetration;
        if (dx * this.hor > 0 && dy * this.ver > 0) {
            const dist = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
            x = this.xpos + this.radius * dx / dist;
            y = this.ypos + this.radius * dy / dist;
            if (this.convex)
                penetration = radius + this.radius - dist;
            else
                penetration = dist + radius - this.radius;
        } else {
            const [hor, ver] = this.convex ? [-this.hor, -this.ver] : [this.hor, this.ver];
            if (dx * hor <= 0 && dy * ver > 0) {
                x = this.xpos;
                y = this.ypos + this.radius * ver;
            } else if (dx * hor > 0 && dy * ver <= 0) {
                x = this.xpos + this.radius * hor;
                y = this.ypos;
            } else return [[this.xpos, this.ypos], 0];
            penetration = 10 - Math.sqrt(Math.pow(xpos - x, 2) + Math.pow(ypos - y, 2));
        }
        if (penetration < 0.0000001 || penetration > 10)
            penetration = 0;
        return [[x, y], penetration];
    }

    // Return the time of intersection (as a fraction of a frame) for the closest point in the ninja's path.
    // Return 0 if the ninja is already intersecting or 1 if it won't intersect within the frame.
    intersect_with_ray(xpos, ypos, dx, dy, radius) {
        // TODO
    }
}

// Define a multidimensional array
// NOTE: Do not try shortening this function by removing the map() method, or else bad things will happen
const xy_array = (x, y, value) => new Array(x).fill(0).map(() => new Array(y).fill(value));

// Return the cartesian product of two number ranges
function product(ranges) {
    for (let i = 0; i < ranges.length; i++) {
        const [start, end] = [ranges[i][0], ranges[i][1]];
        ranges[i] = new Array(end - start).fill(start).map((v, i) => v + i);
    }
    return ranges[0].map(x => ranges[1].map(y => [x, y])).flat();
}

// Return a list containing all cells that contain a rectangle bounded by 2 points
function gather_cells_from_region(x1, y1, x2, y2) {
    const cx1 = Math.floor(x1 / 24);
    const cy1 = Math.floor(y1 / 24);
    const cx2 = Math.floor(x2 / 24);
    const cy2 = Math.floor(y2 / 24);
    return product([[cx1, cx2 + 1], [cy1, cy2 + 1]]);
}

// Return a list that contains all the entities from the nine neighbour cells
function gather_neighbour_cells_content(xpos, ypos) {
    const cx = Math.floor(xpos / 24);
    const cy = Math.floor(ypos / 24);
    const cells = product([[cx - 1, cx + 2], [cy - 1, cy + 2]]);
    return cells.map(v => entity_dic[v]);
}

// Fetch all segments from neighbourhood. If there is collision, return the closest intersection point.
function get_single_closest_point(xpos, ypos, radius) {
    const neighbour_cells = gather_cells_from_region(xpos - radius, ypos - radius, xpos + radius, ypos + radius);
    let biggest_penetration = 0;
    for (const cell of neighbour_cells) {
        for (const segment of segment_dic[cell[0]][cell[1]]) {
            if (segment.active) {
                const [point, penetration] = segment.collision_check(xpos, ypos, radius);
                if (penetration > biggest_penetration) {
                    biggest_penetration = penetration;
                    closest_point = point;
                }
            }
        }
    }
    if (biggest_penetration > 0) return closest_point;
}

// Fetch all segments from neighbourhood. Return shortest intersection time from interpolation.
function sweep_circle_vs_tiles(xpos_old, ypos_old, dx, dy, radius) {
    const xpos_new = xpos_old + dx;
    const ypos_new = ypos_old + dy;
    const width = radius + 1;
    const x1 = Math.min(xpos_old, xpos_new) - width;
    const y1 = Math.min(ypos_old, ypos_new) - width;
    const x2 = Math.max(xpos_old, xpos_new) + width;
    const y2 = Math.max(ypos_old, ypos_new) + width;
    const cells = gather_cells_from_region(x1, y1, x2, y2);
    let shortest_time = 1;
    for (const cell of cells) {
        for (const segment of segment_dic[cell[0]][cell[1]]) {
            if (segment.active && segment.type == "linear") {
                const time = segment.intersect_with_ray(xpos_old, ypos_old, dx, dy, radius);
                shortest_time = Math.min(time, shortest_time);
            }
        }
    }
    return shortest_time;
}

// Return time of intersection by interpolation by sweeping a circle onto an other circle, given a combined radius
function get_time_of_intersection_circle_vs_circle(xpos, ypos, vx, vy, a, b, radius) {
    const dx = xpos - a;
    const dy = ypos - b;
    const dist_sq = Math.pow(dx, 2) + Math.pow(dy, 2);
    const radius_sq = Math.pow(radius, 2);
    const vel_sq = Math.pow(vx, 2) + Math.pow(vy, 2);
    const dot_prod = dx * vx + dy * vy;
    if (dist_sq - radius_sq > 0) {
        if (vel_sq > 0.0001 && dot_prod < 0 && Math.pow(dot_prod, 2) >= vel_sq * (dist_sq - radius_sq))
            return (-dot_prod - Math.sqrt(Math.pow(dot_prod, 2) - vel_sq * (dist_sq - radius_sq))) / vel_sq;
        return 1;
    }
    return 0;
}

// Return time of intersection by interpolation by sweeping a circle onto a line segment
function get_time_of_intersection_circle_vs_lineseg(xpos, ypos, dx, dy, a1, b1, a2, b2, radius) {
    const wx = a2 - a1;
    const wy = b2 - b1;
    const seg_len = Math.sqrt(Math.pow(wx, 2) + Math.pow(wy, 2));
    const nx = wx / seg_len;
    const ny = wy / seg_len;
    const normal_proj = (xpos - a1) * ny - (ypos - b1) * nx;
    const hor_proj = (xpos - a1) * nx + (ypos - b1) * ny;
    if (Math.abs(normal_proj) >= radius) {
        const dir = dx * ny - dy * nx;
        if (dir * normal_proj < 0) {
            const t = Math.min(Math.abs(normal_proj - radius) / Math.abs(dir), 1);
            const hor_proj2 = hor_proj + t * (dx * nx + dy * ny);
            if (0 <= hor_proj2 && hor_proj2 <= seg_len) return t;
        }
    } else if (0 <= hor_proj && hor_proj <= seg_len) return 0;
    return 1;
}

// Return time of intersection by interpolation by sweeping a circle onto a circle arc
function get_time_of_intersection_circle_vs_arc() {
    // TODO
}

// Return the depenetration vector to depenetrate a point out of a square.
// Used to collide the ninja with square entities. (bounce blocks, thwumps, shwumps)
function collision_square_vs_point(square_pos, point_pos, semiside, radius) {
    const x0 = square_pos[0];
    const y0 = square_pos[1];
    const x1 = point_pos[0];
    const y1 = point_pos[1];
    const dx = x1 - x0;
    const dy = y1 - y0;
    const penx = semiside + radius - Math.abs(dx);
    const peny = semiside + radius - Math.abs(dy);
    if (penx > 0 && peny > 0) {
        if (peny <= penx)
            return dy < 0 ? [0, -peny] : [0, peny];
        return dx < 0 ? [-penx, 0] : [penx, 0];
    }
    return [0, 0];
}

// Return a normalized vector pointing in the direction of the orientation.
// Orientation is a value between 0 and 7 taken from map data.
function map_orientation_to_vector(orientation) {
    const diag = Math.sqrt(2) / 2;
    const orientation_dic = [[1, 0], [diag, diag], [0, 1], [-diag, diag], [-1, 0], [-diag, -diag], [0, -1], [diag, -diag]];
    return orientation_dic[orientation];
}

// Return true if the cell has no solid horizontal edge in the specified direction
function is_empty_row(xcoord1, xcoord2, ycoord, dir) {
    const xcoord3 = xcoord1 == xcoord2 ? xcoord1 : xcoord1 + 1;
    if (dir == 1)
        return !(hor_grid_edge_dic[xcoord1][ycoord + 1] || hor_grid_edge_dic[xcoord2][ycoord + 1] || hor_grid_edge_dic[xcoord3][ycoord + 1]);
    if (dir == -1)
        return !(hor_grid_edge_dic[xcoord1][ycoord] || hor_grid_edge_dic[xcoord2][ycoord] || hor_grid_edge_dic[xcoord3][ycoord]);
}

// Return true if the cell has no solid vertical edge in the specified direction
function is_empty_column(xcoord, ycoord1, ycoord2, dir) {
    const ycoord3 = ycoord1 == ycoord2 ? ycoord1 : ycoord1 + 1;
    if (dir == 1)
        return !(ver_grid_edge_dic[xcoord + 1][ycoord1] || ver_grid_edge_dic[xcoord + 1][ycoord2] || ver_grid_edge_dic[xcoord + 1][ycoord3]);
    if (dir == -1)
        return !(ver_grid_edge_dic[xcoord][ycoord1] || ver_grid_edge_dic[xcoord][ycoord2] || ver_grid_edge_dic[xcoord][ycoord3]);
}

// Draw to canvas
function draw() {
    // Prepare canvas for current frame
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Draw tiles
    ctx.fillStyle = "#797988";
    ctx.textAlign = "center";
    for (let x = 0; x < tile_dic.length; x++) {
        for (let y = 0; y < tile_dic[x].length; y++) {
            const tile = tile_dic[x][y];
            if (tile == 1 || tile > 33) {
                ctx.fillRect(x * 24, y * 24, 24, 24);
            } else if (tile > 1) {
                if (tile < 6) {
                    const dx = tile == 3 ? 12 : 0;
                    const dy = tile == 4 ? 12 : 0;
                    const w = tile % 2 == 0 ? 24 : 12;
                    const h = tile % 2 == 0 ? 12 : 24;
                    ctx.fillRect(x * 24 + dx, y * 24 + dy, w, h);
                } else if (tile < 10) {
                    const dx1 = 0;
                    const dy1 = tile == 8 ? 24 : 0;
                    const dx2 = tile == 9 ? 0 : 24;
                    const dy2 = tile == 9 ? 24 : 0;
                    const dx3 = tile == 6 ? 0 : 24;
                    const dy3 = 24;
                    ctx.beginPath();
                    ctx.moveTo(x * 24 + dx1, y * 24 + dy1);
                    ctx.lineTo(x * 24 + dx2, y * 24 + dy2);
                    ctx.lineTo(x * 24 + dx3, y * 24 + dy3);
                    ctx.closePath();
                    ctx.fill();
                } else if (tile < 14) {
                    const dx = (tile == 11 || tile == 12) ? 24 : 0;
                    const dy = (tile == 12 || tile == 13) ? 24 : 0;
                    const a1 = (Math.PI / 2) * (tile - 10);
                    const a2 = (Math.PI / 2) * (tile - 9);
                    ctx.beginPath();
                    ctx.moveTo(x * 24 + dx, y * 24 + dy);
                    ctx.arc(x * 24 + dx, y * 24 + dy, 24, a1, a2);
                    ctx.lineTo(x * 24 + dx, y * 24 + dy);
                    ctx.closePath();
                    ctx.fill();
                } else if (tile < 18) {
                    const dx1 = (tile == 15 || tile == 16) ? 24 : 0;
                    const dy1 = (tile == 16 || tile == 17) ? 24 : 0;
                    const dx2 = (tile == 14 || tile == 17) ? 24 : 0;
                    const dy2 = (tile == 14 || tile == 15) ? 24 : 0;
                    const a1 = Math.PI + (Math.PI / 2) * (tile - 10);
                    const a2 = Math.PI + (Math.PI / 2) * (tile - 9);
                    ctx.beginPath();
                    ctx.moveTo(x * 24 + dx1, y * 24 + dy1);
                    ctx.arc(x * 24 + dx2, y * 24 + dy2, 24, a1, a2);
                    ctx.lineTo(x * 24 + dx1, y * 24 + dy1);
                    ctx.closePath();
                    ctx.fill();
                } else if (tile < 22) {
                    const dx1 = 0;
                    const dy1 = (tile == 20 || tile == 21) ? 24 : 0;
                    const dx2 = 24
                    const dy2 = (tile == 20 || tile == 21) ? 24 : 0;
                    const dx3 = (tile == 19 || tile == 20) ? 24 : 0;
                    const dy3 = 12;
                    ctx.beginPath();
                    ctx.moveTo(x * 24 + dx1, y * 24 + dy1);
                    ctx.lineTo(x * 24 + dx2, y * 24 + dy2);
                    ctx.lineTo(x * 24 + dx3, y * 24 + dy3);
                    ctx.closePath();
                    ctx.fill();
                } else if (tile < 26) {
                    const dx1 = 0;
                    const dy1 = (tile == 23 || tile == 24) ? 12 : 0;
                    const dx2 = tile == 23 ? 0 : 24;
                    const dy2 = tile == 25 ? 12 : 0;
                    const dx3 = 24;
                    const dy3 = tile < 24 ? (tile == 22 ? 12 : 0) : 24;
                    const dx4 = tile == 23 ? 24 : 0;
                    const dy4 = 24;
                    ctx.beginPath();
                    ctx.moveTo(x * 24 + dx1, y * 24 + dy1);
                    ctx.lineTo(x * 24 + dx2, y * 24 + dy2);
                    ctx.lineTo(x * 24 + dx3, y * 24 + dy3);
                    ctx.lineTo(x * 24 + dx4, y * 24 + dy4);
                    ctx.closePath();
                    ctx.fill();
                } else if (tile < 30) {
                    const dx1 = 12;
                    const dy1 = (tile == 28 || tile == 29) ? 24 : 0;
                    const dx2 = (tile == 27 || tile == 28) ? 24 : 0
                    const dy2 = 0;
                    const dx3 = (tile == 27 || tile == 28) ? 24 : 0;
                    const dy3 = 24;
                    ctx.beginPath();
                    ctx.moveTo(x * 24 + dx1, y * 24 + dy1);
                    ctx.lineTo(x * 24 + dx2, y * 24 + dy2);
                    ctx.lineTo(x * 24 + dx3, y * 24 + dy3);
                    ctx.closePath();
                    ctx.fill();
                } else if (tile < 34) {
                    const dx1 = 12;
                    const dy1 = (tile == 30 || tile == 31) ? 24 : 0;
                    const dx2 = (tile == 31 || tile == 33) ? 24 : 0;
                    const dy2 = 24
                    const dx3 = (tile == 31 || tile == 32) ? 24 : 0;
                    const dy3 = (tile == 32 || tile == 33) ? 24 : 0;
                    const dx4 = (tile == 30 || tile == 32) ? 24 : 0;
                    const dy4 = 0;
                    ctx.beginPath();
                    ctx.moveTo(x * 24 + dx1, y * 24 + dy1);
                    ctx.lineTo(x * 24 + dx2, y * 24 + dy2);
                    ctx.lineTo(x * 24 + dx3, y * 24 + dy3);
                    ctx.lineTo(x * 24 + dx4, y * 24 + dy4);
                    ctx.closePath();
                    ctx.fill();
                }
            }
        }
    }

    // Draw ninja
    ctx.fillStyle = "#000";
    ctx.beginPath();
    ctx.arc(p.xpos, p.ypos, p.radius, 0, Math.PI * 2);
    ctx.fill();
    ctx.closePath();
}

// This is the main function that handles physics.
function tick() {
    // Only perform operations on interval defined by framerate
    const current_frame = Date.now();
    if (unlimited || current_frame - last_frame >= interval) {
        last_frame = current_frame;
        
        // Extract inputs for this frame
        p.hor_input = -input.left + input.right;
        p.jump_input = input.jump;
    
        // TODO
        
        p.integrate(); // Do preliminary speed and position updates
        p.pre_collision(); // Do pre-collision calculations
        for (let i = 0; i < 4; i++) {
            p.collide_vs_objects(); // Handle physical collisions with entities
            p.collide_vs_tiles(); // Handle physical collisions with tiles
        }
        p.post_collision(); // Do post-collision calculations
        p.think(); // Make ninja think
    
        draw(); // Update canvas
    }
    requestAnimationFrame(tick);
}

let initialized = false;
function init(mdata_buffer) {
    // Import map data
    mdata = new Int8Array(mdata_buffer);

    // Extract tile data from map data
    tile_data = mdata.slice(184, 1150);

    // Initiate a dictionary mapping each tile to its cell. Start by filling it with full tiles.
    tile_dic = xy_array(44, 25, 1);

    // Map each tile to its cell, leaving the border tiles solid
    for (let x = 0; x < tile_dic.length - 2; x++)
        for (let y = 0; y < tile_dic[x].length - 2; y++)
            tile_dic[x + 1][y + 1] = tile_data[x + y * 42];
    
    // Initiate dictionaries and list containing interactable segments and entities
    segment_dic = xy_array(45, 26, []);
    entity_dic = xy_array(45, 26, []);
    entity_list = [];

    // Initiate dictionaries of grid edges and segments. They are all set to zero initially.
    // Increment by 1 whenever a grid edge or segment is spawned.
    hor_grid_edge_dic = xy_array(88, 51, 0);
    ver_grid_edge_dic = xy_array(89, 50, 0);
    hor_segment_dic = xy_array(88, 51, 0);
    ver_segment_dic = xy_array(89, 50, 0);

    tile_grid_edge_map = [
        [0,0,0,0,0,0,0,0,0,0,0,0], [1,1,0,0,1,1,1,1,0,0,1,1], // Empty and full tiles
        [1,1,1,1,0,0,1,0,0,0,1,0], [0,1,0,0,0,1,0,0,1,1,1,1], [0,0,1,1,1,1,0,1,0,0,0,1], [1,0,0,0,1,0,1,1,1,1,0,0], // Half tiles
        [1,1,0,1,1,0,1,1,0,1,1,0], [1,1,1,0,0,1,1,0,0,1,1,1], [0,1,1,0,1,1,0,1,1,0,1,1], [1,0,0,1,1,1,1,1,1,0,0,1], // 45 degree slopes
        [1,1,0,0,1,1,1,1,0,0,1,1], [1,1,0,0,1,1,1,1,0,0,1,1], [1,1,0,0,1,1,1,1,0,0,1,1], [1,1,0,0,1,1,1,1,0,0,1,1], // Quarter moons
        [1,1,0,1,1,0,1,1,0,1,1,0], [1,1,1,0,0,1,1,0,0,1,1,1], [0,1,1,0,1,1,0,1,1,0,1,1], [1,0,0,1,1,1,1,1,1,0,0,1], // Quarter pipes
        [1,1,1,1,0,0,1,0,0,0,1,0], [1,1,1,1,0,0,1,0,0,0,1,0], [0,0,1,1,1,1,0,1,0,0,0,1], [0,0,1,1,1,1,0,1,0,0,0,1], // Short mild slopes
        [1,1,0,0,1,1,1,1,0,0,1,1], [1,1,0,0,1,1,1,1,0,0,1,1], [1,1,0,0,1,1,1,1,0,0,1,1], [1,1,0,0,1,1,1,1,0,0,1,1], // Raised mild slopes
        [1,0,0,0,1,0,1,1,1,1,0,0], [0,1,0,0,0,1,0,0,1,1,1,1], [0,1,0,0,0,1,0,0,1,1,1,1], [1,0,0,0,1,0,1,1,1,1,0,0], // Short steep slopes
        [1,1,0,0,1,1,1,1,0,0,1,1], [1,1,0,0,1,1,1,1,0,0,1,1], [1,1,0,0,1,1,1,1,0,0,1,1], [1,1,0,0,1,1,1,1,0,0,1,1], // Raised steep slopes
        [1,1,0,0,0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0,0,1,1], [0,0,0,0,1,1,0,0,0,0,0,0], [0,0,0,0,0,0,1,1,0,0,0,0]  // Glitched tiles
    ];

    tile_segment_map = [
        [0,0,0,0,0,0,0,0,0,0,0,0], [1,1,0,0,1,1,1,1,0,0,1,1], // Empty and full tiles
        [1,1,1,1,0,0,1,0,0,0,1,0], [0,1,0,0,0,1,0,0,1,1,1,1], [0,0,1,1,1,1,0,1,0,0,0,1], [1,0,0,0,1,0,1,1,1,1,0,0], // Half tiles
        [1,1,0,0,0,0,1,1,0,0,0,0], [1,1,0,0,0,0,0,0,0,0,1,1], [0,0,0,0,1,1,0,0,0,0,1,1], [0,0,0,0,1,1,1,1,0,0,0,0], // 45 degree slopes
        [1,1,0,0,0,0,1,1,0,0,0,0], [1,1,0,0,0,0,0,0,0,0,1,1], [0,0,0,0,1,1,0,0,0,0,1,1], [0,0,0,0,1,1,1,1,0,0,0,0], // Quarter moons
        [1,1,0,0,0,0,1,1,0,0,0,0], [1,1,0,0,0,0,0,0,0,0,1,1], [0,0,0,0,1,1,0,0,0,0,1,1], [0,0,0,0,1,1,1,1,0,0,0,0], // Quarter pipes
        [1,1,0,0,0,0,1,0,0,0,0,0], [1,1,0,0,0,0,0,0,0,0,1,0], [0,0,0,0,1,1,0,0,0,0,0,1], [0,0,0,0,1,1,0,1,0,0,0,0], // Short mild slopes
        [1,1,0,0,0,0,1,1,0,0,1,0], [1,1,0,0,0,0,1,0,0,0,1,1], [0,0,0,0,1,1,0,1,0,0,1,1], [0,0,0,0,1,1,1,1,0,0,0,1], // Raised mild slopes
        [1,0,0,0,0,0,1,1,0,0,0,0], [0,1,0,0,0,0,0,0,0,0,1,1], [0,0,0,0,0,1,0,0,0,0,1,1], [0,0,0,0,1,0,1,1,0,0,0,0], // Short steep slopes
        [1,1,0,0,1,0,1,1,0,0,0,0], [1,1,0,0,0,1,0,0,0,0,1,1], [0,1,0,0,1,1,0,0,0,0,1,1], [1,0,0,0,1,1,1,1,0,0,0,0], // Raised steep slopes
        [1,1,0,0,0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0,0,1,1], [0,0,0,0,1,1,0,0,0,0,0,0], [0,0,0,0,0,0,1,1,0,0,0,0]  // Glitched tiles
    ];
    
    for (let xcoord = 0; xcoord < tile_dic.length; xcoord++) {
        for (let ycoord = 0; ycoord < tile_dic[xcoord].length; ycoord++) {
            let tile_id = tile_dic[xcoord][ycoord];
            if (tile_id > 37) tile_id = 0;
            const grid_edge_list = tile_grid_edge_map[tile_id];
            const segment_list = tile_segment_map[tile_id];
            for (let y = 0; y < 3; y++) {
                for (let x = 0; x < 2; x++) {
                    hor_grid_edge_dic[2 * xcoord + x][2 * ycoord + y] += grid_edge_list[2 * y + x];
                    hor_segment_dic[2 * xcoord + x][2 * ycoord + y] += segment_list[2 * y + x];
                }
            }
            for (let x = 0; x < 3; x++) {
                for (let y = 0; y < 2; y++) {
                    ver_grid_edge_dic[2 * xcoord + x][2 * ycoord + y] += grid_edge_list[2 * x + y + 6];
                    ver_segment_dic[2 * xcoord + x][2 * ycoord + y] += segment_list[2 * x + y + 6];
                }
            }
            const xtl = xcoord * 24;
            const ytl = ycoord * 24;
            if (tile_id == 6 || tile_id == 8)
                segment_dic[xcoord][ycoord].push(new GridSegmentLinear([xtl, ytl + 24], [xtl + 24, ytl]));
            else if (tile_id == 7 || tile_id == 9)
                segment_dic[xcoord][ycoord].push(new GridSegmentLinear([xtl, ytl], [xtl + 24, ytl + 24]));
            else if (tile_id == 10)
                segment_dic[xcoord][ycoord].push(new GridSegmentCircular([xtl, ytl], [1, 1]));
            else if (tile_id == 11)
                segment_dic[xcoord][ycoord].push(new GridSegmentCircular([xtl + 24, ytl], [-1, 1]));
            else if (tile_id == 12)
                segment_dic[xcoord][ycoord].push(new GridSegmentCircular([xtl + 24, ytl + 24], [-1, -1]));
            else if (tile_id == 13)
                segment_dic[xcoord][ycoord].push(new GridSegmentCircular([xtl, ytl + 24], [1, -1]));
            else if (tile_id == 14)
                segment_dic[xcoord][ycoord].push(new GridSegmentCircular([xtl + 24, ytl + 24], [-1, -1], false));
            else if (tile_id == 15)
                segment_dic[xcoord][ycoord].push(new GridSegmentCircular([xtl, ytl + 24], [1, -1], false));
            else if (tile_id == 16)
                segment_dic[xcoord][ycoord].push(new GridSegmentCircular([xtl, ytl], [1, 1], false));
            else if (tile_id == 17)
                segment_dic[xcoord][ycoord].push(new GridSegmentCircular([xtl + 24, ytl], [-1, 1], false));
            else if (tile_id == 18 || tile_id == 24)
                segment_dic[xcoord][ycoord].push(new GridSegmentLinear([xtl, ytl + 12], [xtl + 24, ytl]));
            else if (tile_id == 19 || tile_id == 25)
                segment_dic[xcoord][ycoord].push(new GridSegmentLinear([xtl, ytl], [xtl + 24, ytl + 12]));
            else if (tile_id == 20 || tile_id == 22)
                segment_dic[xcoord][ycoord].push(new GridSegmentLinear([xtl, ytl + 24], [xtl + 24, ytl + 12]));
            else if (tile_id == 21 || tile_id == 23)
                segment_dic[xcoord][ycoord].push(new GridSegmentLinear([xtl, ytl + 12], [xtl + 24, ytl + 24]));
            else if (tile_id == 26 || tile_id == 32)
                segment_dic[xcoord][ycoord].push(new GridSegmentLinear([xtl, ytl + 24], [xtl + 12, ytl]));
            else if (tile_id == 27 || tile_id == 33)
                segment_dic[xcoord][ycoord].push(new GridSegmentLinear([xtl + 12, ytl], [xtl + 24, ytl + 24]));
            else if (tile_id == 28 || tile_id == 30)
                segment_dic[xcoord][ycoord].push(new GridSegmentLinear([xtl + 12, ytl + 24], [xtl + 24, ytl]));
            else if (tile_id == 29 || tile_id == 31)
                segment_dic[xcoord][ycoord].push(new GridSegmentLinear([xtl, ytl], [xtl + 12, ytl + 24]));
        }
    }
    for (let xcoord = 0; xcoord < hor_segment_dic.length; xcoord++) {
        for (let ycoord = 0; ycoord < hor_segment_dic[xcoord].length; ycoord++) {
            const state = hor_segment_dic[xcoord][ycoord];
            if (state == 1) segment_dic[Math.floor(xcoord / 2)][Math.floor(ycoord / 2)].push(new GridSegmentLinear([12 * xcoord, 12 * ycoord], [12 * xcoord + 12, 12 * ycoord]));
        }
    }
    for (let xcoord = 0; xcoord < ver_segment_dic.length; xcoord++) {
        for (let ycoord = 0; ycoord < ver_segment_dic[xcoord].length; ycoord++) {
            const state = ver_segment_dic[xcoord][ycoord];
            if (state == 1) segment_dic[Math.floor(xcoord / 2)][Math.floor(ycoord / 2)].push(new GridSegmentLinear([12 * xcoord, 12 * ycoord], [12 * xcoord, 12 * ycoord + 12]));
        }
    }
    
    // Find the spawn position of the ninja
    const xspawn = mdata[1231] * 6;
    const yspawn = mdata[1232] * 6;

    // Initiate player 1 instance of ninja at spawn coordinates
    p = new Ninja(xspawn, yspawn);

    // TODO

    // Execute the main physics function once per frame
    if (!initialized) {
        initialized = true;
        tick();
    }
}

// Load placeholder map upon page load
(async () => init(await fetch("default_map").then(r => r.blob()).then(r => r.arrayBuffer())))();

// Load map from computer
document.querySelector("input[type=file]").addEventListener("change", e => {
    const file = e.target.files[0];
    const reader = new FileReader();
    reader.addEventListener("load", e => init(e.target.result));
    reader.readAsArrayBuffer(file);
});

// Modify framerate
document.querySelector("input[type=number]").value = fps;
document.querySelector("input[type=number]").addEventListener("change", e => {
    const new_fps = parseInt(e.target.value);
    if (isNaN(new_fps)) return;
    fps = new_fps;
    interval = 1000 / fps;
});

// Toggle unlimited framerate
document.querySelector("input[type=checkbox]").addEventListener("change", e => {
    unlimited = e.target.checked;
});