import {orient2d} from 'robust-predicates';
import containingTriangle from '@kninnug/containing-triangle';

export interface DelaunatorLike {
    coords: {readonly length: number, readonly [n: number]: number};
    triangles:  {readonly length: number, [n: number]: number};
    halfedges:  {readonly length: number, [n: number]: number};
    hull:  {readonly length: number, readonly [n: number]: number};
};

/** A 4-tuple indicating the left-x/-y and right-x/-y of a segment. */
export type Segment = [number, number, number, number];

/**
 * Compute if and where a ray, or half-line, intersects with a segment. Returns
 * the relative location of the intersection point on the ray. I.e, if it
 * intersects:
 *
 *     const r = segIntersectRay(s1x, s1y, s2x, s2y, r1x, r1y, r2x, r2y),
 *           dx = r2x - r1x,
 *           dy = r2y - r1y,
 *           // intersection is at:
 *           x = r1x + r * dx,
 *           y = r1y + r * dy;
 *
 * @source https://ncase.me/sight-and-light/
 * @param s1x The x coordinate of point 1 of the segment.
 * @param s1y The y coordinate of point 1 of the segment.
 * @param s2x The x coordinate of point 2 of the segment.
 * @param s2y The y coordinate of point 2 of the segment.
 * @param r1x The x coordinate of the origin of the ray.
 * @param r1y The y coordinate of the origin of the ray.
 * @param r2x The x coordinate of a second point along the ray.
 * @param r2y The y coordinate of a second point along the ray.
 * @return The relative point of the intersection, or Infinity if the
 *         ray never intersects the segment.
 */
function segIntersectRay(s1x: number, s1y: number, s2x: number, s2y: number, r1x: number, r1y: number, r2x: number, r2y: number){
    const rdx = r2x - r1x,
        rdy = r2y - r1y,
        
        sdx = s2x - s1x,
        sdy = s2y - s1y,
        
        rmag = Math.sqrt(rdx * rdx + rdy * rdy),
        smag = Math.sqrt(sdx * sdx + sdy * sdy);
    
    if(rdx / rmag === sdx / smag && rdy / rmag === sdy / smag){
        return Infinity;
    }
    
    const T2 = (rdx * (s1y - r1y) + rdy * (r1x - s1x)) / (sdx * rdy - sdy * rdx),
        T1 = rdx ? (s1x + sdx*T2 - r1x) / rdx : (s1y + sdy*T2 - r1y) / rdy;
    
    if(T1 < 0 || T2 < 0 || T2 > 1.0){
        return Infinity;
    }
    
    return T1;
}

/**
 * Square distance between two points of a triangulation.
 *
 * @param x1 First point x coordinate.
 * @param y1 First point y coordinate.
 * @param x2 Second point x coordinate.
 * @param y2 Second point y coordinate.
 * @return The squared distance.
 */
function sqdist(x1: number, y1: number, x2: number, y2: number){
    const dx = x2 - x1,
        dy = y2 - y1;
    
    return dx * dx + dy * dy;
}

/**
 * Next half-edge counter-clockwise in a triangle.
 *
 * @param e Half-edge id.
 * @return Id of the next half-edge.
 */
function nextEdge(e: number){ return (e % 3 === 2) ? e - 2 : e + 1; }
/**
 * Previous half-edge counter-clockwise in a triangle.
 *
 * @param e Half-edge id.
 * @return Id of the previous half-edge.
 */
function prevEdge(e: number){ return (e % 3 === 0) ? e + 2 : e - 1; }
/**
 * Half-edges of the given triangle.
 *
 * @param t Triangle id.
 * @return Ids of the half-edges.
 */
function edgesOfTri(t: number){ return [t * 3, t * 3 + 1, t * 3 + 2]; }
/**
 * Point indices of the given triangle.
 *
 * @param del The triangulation.
 * @param t Triangle id.
 * @return Indices into the points array of the triangle's points.
 */
function pointsOfTri(del: DelaunatorLike, t: number){
    const p1 = del.triangles[t * 3],
        p2 = del.triangles[t * 3 + 1],
        p3 = del.triangles[t * 3 + 2];
    return [p1, p2, p3];
}

/**
 * Whether a point is left of the line defined by two other points.
 *
 * @param x1 The x coordinate of the first point on the line.
 * @param y1 The y coordinate of the first point on the line.
 * @param x2 The x coordinate of the second point on the line.
 * @param y2 The y coordinate of the second point on the line.
 * @param px The x coordinate of the query point.
 * @param py The y coordinate of the query point.
 * @return True if (px, py) is strictly to the left of the segment.
 */
function isLeftOf(x1: number, y1: number, x2: number, y2: number, px: number, py: number){
    return orient2d(x1, y1, x2, y2, px, py) > 0;
}

/**
 * Whether a point is right of the line defined by two other points.
 *
 * @param x1 The x coordinate of the first point on the line.
 * @param y1 The y coordinate of the first point on the line.
 * @param x2 The x coordinate of the second point on the line.
 * @param y2 The y coordinate of the second point on the line.
 * @param px The x coordinate of the query point.
 * @param py The y coordinate of the query point.
 * @return True if (px, py) is strictly to the right of the segment.
 */
function isRightOf(x1: number, y1: number, x2: number, y2: number, px: number, py: number){
    return orient2d(x1, y1, x2, y2, px, py) < 0;
}

/**
 * Compute the centroid of a triangle in a triangulation.
 *
 * @param del The triangulation.
 * @param tri The triangle id.
 * @return The coordinates of the centroid.
 */
function centroid(del: DelaunatorLike, tri: number){
    const [p1, p2, p3] = pointsOfTri(del, tri),
        p1x = del.coords[p1 * 2],
        p1y = del.coords[p1 * 2 + 1],
        p2x = del.coords[p2 * 2],
        p2y = del.coords[p2 * 2 + 1],
        p3x = del.coords[p3 * 2],
        p3y = del.coords[p3 * 2 + 1];
    return [(p1x + p2x + p3x) / 3, (p1y + p2y + p3y) / 3];
}

/**
 * Order two points with respect to a query point.
 *
 * @param qx The x coordinate of the query point.
 * @param qy The y coordinate of the query point.
 * @param p1x The x coordinate of the first point.
 * @param p1y The y coordinate of the first point.
 * @param p2x The x coordinate of the second point.
 * @param p2y The y coordinate of the second point.
 * @return The coordinates of the two points in order [left-x,
 *         left-y, right-x, right-y].
 */
function orderAngles(qx: number, qy: number, p1x: number, p1y: number, p2x: number, p2y: number){
    const segLeft = isLeftOf(qx, qy, p2x, p2y, p1x, p1y),
        segRight = !segLeft,
        lx = segLeft ? p1x : p2x,
        ly = segLeft ? p1y : p2y,
        rx = segRight ? p1x : p2x,
        ry = segRight ? p1y : p2y;
    return [lx, ly, rx, ry];
}

/**
 * Order two points of a triangulation with respect to a query point.
 *
 * @param del The triangulation.
 * @param qx The x coordinate of the query point.
 * @param qy The y coordinate of the query point.
 * @param p1 The id of the first point.
 * @param p2 The id of the second point.
 * @return The coordinates of the two points in order [left-x,
 *         left-y, right-x, right-y].
 */
function orderDelAngles(del: DelaunatorLike, qx: number, qy: number, p1: number, p2: number){
    const p1x = del.coords[p1 * 2],
        p1y = del.coords[p1 * 2 + 1],
        p2x = del.coords[p2 * 2],
        p2y = del.coords[p2 * 2 + 1];
    return orderAngles(qx, qy, p1x, p1y, p2x, p2y);
}

/**
 * Whether a segment is within a viewing cone.
 *
 * @param px The x coordinate of the query point.
 * @param py The y coordinate of the query point.
 * @param slx The x coordinate of the segment's left point.
 * @param sly The y coordinate of the segment's left point.
 * @param srx The x coordinate of the segment's right point.
 * @param sry The y coordinate of the segment's right point.
 * @param rlx The x coordinate of the viewing cone's left point.
 * @param rly The y coordinate of the viewing cone's left point.
 * @param rrx The x coordinate of the viewing cone's right point.
 * @param rry The y coordinate of the viewing cone's right point.
 * @return True if the segment is within the viewing cone.
 */
function isWithinCone(px: number, py: number, slx: number, sly: number, srx: number, sry: number, rlx: number, rly: number, rrx: number, rry: number){
    //  o >----- o #====# o -----< o
    //  rl       sl      sr       rr  -->  (rl leftOf sl && sl leftOf rr) || (rl leftOf sr && sr leftOf rr)
    //  o #===== o >====# o -----< o
    //  sl       rl      sr       rr  -->  sl leftOr rl && rl leftOf sr
    //  o >----- o #====< o =====# o
    //  rl       sl      rr       sr  -->  rl leftOf sl && sl leftOf rr
    //  o #===== o >====< o =====# o
    //  sl       rl      rr       sr  -->  (sl leftOf rl && rl leftOf sr) || (sl leftOf rr && rr leftOf sr)
    
    //  o >----< o ------ o #====# o
    //  rl      rr        sl      sr  -->  rr leftOf sl
    if(isLeftOf(px, py, slx, sly, rrx, rry)){
        return false;
    }
    //  o #====# o ------ o >----< o
    //  sl      sr        rl      rr  -->  sr leftOf rl
    if(isLeftOf(px, py, rlx, rly, srx, sry)){
        return false;
    }
    //  o >--------< o #========# o
    //  rl         rr sl         sr  -->  
    if(rrx === slx && rry === sly){
        return false;
    }
    //  o #========# o >--------< o
    //  sl         sr rl         rr  -->  
    if(srx === rlx && sry === rly){
        return false;
    }
    
    return true;
}

/**
 * Restrict a viewing cone by an edge.
 *
 * @param px The x coordinate of the query point.
 * @param py The y coordinate of the query point.
 * @param slx The x coordinate of the segment's left point.
 * @param sly The y coordinate of the segment's left point.
 * @param srx The x coordinate of the segment's right point.
 * @param sry The y coordinate of the segment's right point.
 * @param rlx The x coordinate of the viewing cone's left point.
 * @param rly The y coordinate of the viewing cone's left point.
 * @param rrx The x coordinate of the viewing cone's right point.
 * @param rry The y coordinate of the viewing cone's right point.
 * @return The coordinates of the new viewing cone, as
 *         well as two booleans indicating whether the left and/or the right
 *         points were changed from the original cone to the segment.
 */
function restrictAngles(px: number, py: number,
        slx: number, sly: number, srx: number, sry: number,
        rlx: number, rly: number, rrx: number, rry: number): [...Segment, boolean, boolean] {
    let nlx = rlx,
        nly = rly,
        resLeft = false;
    if(isRightOf(px, py, rlx, rly, slx, sly)){
        nlx = slx;
        nly = sly;
        resLeft = true;
    }
    
    let nrx = rrx,
        nry = rry,
        resRight = false;
    if(isLeftOf(px, py, rrx, rry, srx, sry)){
        nrx = srx;
        nry = sry;
        resRight = true;
    }
    
    return [nlx, nly, nrx, nry, resLeft, resRight];
}

/**
 * Compute a visibility polygon of line segments that are visible from a given
 * query point in a triangulation.
 *
 * @source "Efficient Computation of Visibility Polygons" - Francisc Bungiu et al.
 * @param del The triangulation.
 * @param qx The x coordinate of the query point.
 * @param qy The y coordinate of the query point.
 * @param obstructs A callback that indicates whether an edge in the
 *        triangulation obstructs visibility. Called with an edge id and
 *        expected to return a truthy value.
 * @param ilx The x coordinate of the left point restricting the viewing cone.
 * @param ily The y coordinate of the left point restricting the viewing cone.
 * @param irx The x coordinate of the right point restricting the viewing cone.
 * @param iry The y coordinate of the right point restricting the viewing cone.
 * @return An array of number quadruples that define the
 *         line segments making up the visibility polygon.
 */
function triangularExpansion(del: DelaunatorLike, qx: number, qy: number, obstructs: (edge: number) => boolean,
        ilx = NaN, ily = NaN, irx = NaN, iry = NaN): Segment[] {
    /**
     * Expand the visibility polygon through the next triangle.
     *
     * @param edgIn The edge in this triangle through which we came in.
     * @param rlx The x coordinate of the left restricting point.
     * @param rly The y coordinate of the left restricting point.
     * @param rrx The x coordinate of the right restricting point.
     * @param rry The y coordinate of the right restricting point.
     * @return The segment quadruples of this subsection of
     *         the visibility polygon.
     */
    function expand(edgIn: number, rlx: number, rly: number, rrx: number, rry: number): Segment[] {
        const ret: Segment[] = [],
            edges = [nextEdge(edgIn), prevEdge(edgIn)];
        
        for(let i = 0; i < 2; i++){
            const edg = edges[i],
                p1 = del.triangles[edg],
                p2 = del.triangles[nextEdge(edg)],
                adjOut = del.halfedges[edg];
            // Segment left and right
            let [slx, sly, srx, sry] = orderDelAngles(del, qx, qy, p1, p2);
            
            // Edge outside viewing cone
            if(!isWithinCone(qx, qy, slx, sly, srx, sry, rlx, rly, rrx, rry)){
                continue;
            }
            
            // Next restricting left & right points
            const [nlx, nly, nrx, nry, resL, resR] = restrictAngles(qx, qy,
                    slx, sly, srx, sry, rlx, rly, rrx, rry);
            
            // Viewing cone degenerated to zero-area triangle
            if(orient2d(qx, qy, nrx, nry, nlx, nly) <= 0.0){
                continue;
            }
            
            if(adjOut !== -1 && !obstructs(edg)){
                // Expand into the next triangle
                ret.push(...expand(adjOut, nlx, nly, nrx, nry));
                continue;
            }
            // Hit a wall, or the hull
            
            // Left restricting point is now on the segment
            if(!resL){
                const int = segIntersectRay(slx, sly, srx, sry, qx, qy, rlx, rly);
                if(int !== Infinity){
                    const dx = rlx - qx,
                        dy = rly - qy;
                    slx = qx + int * dx;
                    sly = qy + int * dy;
                }
            }
            // Right restricting point is now on the segment
            if(!resR){
                const int = segIntersectRay(slx, sly, srx, sry, qx, qy, rrx, rry);
                if(int !== Infinity){
                    const dx = rrx - qx,
                        dy = rry - qy;
                    srx = qx + int * dx;
                    sry = qy + int * dy;
                }
            }
            
            ret.push([slx, sly, srx, sry]);
        }
        
        return ret;
    }
    
    const triStart = containingTriangle(del, qx, qy);
    if(triStart === -1){
        return [];
    }
    
    // Order w.r.t. the centroid, rather than the view point to avoid
    // instability when the view point is on a triangle edge.
    const [cx, cy] = centroid(del, triStart),
        startEdges = edgesOfTri(triStart),
        ret: Segment[] = [],
        prestrict = !isNaN(ilx);
    
    // Have an initial viewing cone
    if(prestrict){ // assume all are given, or none
        ([ilx, ily, irx, iry] = orderAngles(qx, qy, ilx, ily, irx, iry));
    }
    
    for(let i = 0; i < 3; i++){
        const edg = startEdges[i],
            p1 = del.triangles[edg],
            p2 = del.triangles[nextEdge(edg)];
        let [slx, sly, srx, sry] = orderDelAngles(del, cx, cy, p1, p2);
        // Fix for stack-overflows when (qx, qy) is exactly on an end-point
        if(sqdist(slx, sly, qx, qy) === 0.0 || sqdist(srx, sry, qx, qy) === 0.0){
            return [];
        }
        
        if(prestrict){
            const intL = segIntersectRay(slx, sly, srx, sry, qx, qy, ilx, ily);
            if(intL !== Infinity){
                const dx = ilx - qx,
                    dy = ily - qy;
                ([, , slx, sly] = orderAngles(qx, qy, qx + intL * dx, qy + intL * dy, slx, sly));
            }
            const intR = segIntersectRay(slx, sly, srx, sry, qx, qy, irx, iry);
            if(intR !== Infinity){
                const dx = irx - qx,
                    dy = iry - qy;
                ([srx, sry, , ] = orderAngles(qx, qy, srx, sry, qx + intR * dx, qy + intR * dy));
            }
            
            ([slx, sly, srx, sry] = orderAngles(qx, qy, slx, sly, srx, sry));
        }else{
            ilx = slx;
            ily = sly;
            irx = srx;
            iry = sry;
        }
        
        const [rlx, rly, rrx, rry] = restrictAngles(qx, qy, slx, sly, srx, sry, ilx, ily, irx, iry);
        
        if(!isWithinCone(qx, qy, slx, sly, srx, sry, ilx, ily, irx, iry)){
            continue;
        }
        
        const adj = del.halfedges[edg];
        if(adj === -1 || obstructs(edg)){
            ret.push([slx, sly, srx, sry]);
        }else{
            ret.push(...expand(adj, rlx, rly, rrx, rry));
        }
    }
    
    return ret;
}

export default triangularExpansion;
