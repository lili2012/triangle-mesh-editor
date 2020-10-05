#include <float.h>
#include <assert.h>
#include "meshEdit.h"
#include "mutablePriorityQueue.h"
#include "error_dialog.h"
#include <set>

struct HalfEdges {
  HalfedgeIter h0;
  HalfedgeIter h1;
};

static bool isInTriangle(HalfedgeIter h)
{
  HalfedgeIter next = h;
  for (int i = 0; i < 3; i++) {
    next = next->next();
  }
  return next == h;
}

static inline int prev(int i, int n)
{
  return (i + n - 1) % n;
}
static inline int next(int i, int n)
{
  return (i + 1) % n;
}

namespace CS248 {
  //see splitEdge.dwg for annotation
  VertexIter HalfedgeMesh::splitEdge(EdgeIter e0) {
    // This method should split the given edge and return an iterator to the
    // newly inserted vertex. The halfedge of this vertex should point along
    // the edge that was split, rather than the new edges.
    // Check e0's face is triangle

    //bool isBoundary = false;
    //halfedge() || halfedge()->twin()->face()->isBoundary();
    if (!checkIfEdgeCouldSplit(e0)) {
      return VertexIter();
    }
    HalfedgeIter h0 = e0->halfedge();
    if (h0->face()->isBoundary()) {
      h0 = h0->twin();
    }

    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h8 = h1->twin();
    HalfedgeIter h9 = h2->twin();

    // VERTICES
    VertexIter v0 = h2->vertex();
    VertexIter v1 = h0->vertex();
    VertexIter v3 = h3->vertex();
    // EDGES
    EdgeIter e3 = h1->edge();
    EdgeIter e4 = h2->edge();
    // FACES
    FaceIter f0 = h0->face();

    //Allocate new elements
    //halfedges
    HalfedgeIter h10 = newHalfedge();
    HalfedgeIter h11 = newHalfedge();
    HalfedgeIter h12 = newHalfedge();
    HalfedgeIter h13 = newHalfedge();

    //vertices
    VertexIter v4 = newVertex();
    v4->position = (v1->position + v3->position) / 2;
    v4->isNew = true;
    //edges
    EdgeIter e5 = newEdge();
    e5->isNew = false;
    EdgeIter e7 = newEdge();
    e7->isNew = true;
    //faces
    FaceIter f3 = newFace();

    //Reassign
    h0->twin() = h10;
    h0->next() = h13;

    h1->next() = h12;
    h1->face() = f3;

    h3->twin() = h11;
    h3->edge() = e5;


    h10->twin() = h0;
    h10->vertex() = v4;
    h10->edge() = e0;

    h11->twin() = h3;
    h11->next() = h1;
    h11->vertex() = v4;
    h11->edge() = e5;
    h11->face() = f3;

    h12->twin() = h13;
    h12->next() = h11;
    h12->vertex() = v0;
    h12->edge() = e7;
    h12->face() = f3;

    h13->twin() = h12;
    h13->next() = h2;
    h13->vertex() = v4;
    h13->edge() = e7;
    h13->face() = f0;

    v4->halfedge() = h10;

    e0->halfedge() = h0;
    e5->halfedge() = h3;
    e7->halfedge() = h12;

    f0->halfedge() = h0;
    f3->halfedge() = h1;

    if (!h3->face()->isBoundary()) {
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();
      HalfedgeIter h6 = h4->twin();
      HalfedgeIter h7 = h5->twin();
      VertexIter v2 = h6->vertex();
      EdgeIter e1 = h4->edge();
      EdgeIter e2 = h5->edge();
      FaceIter f1 = h3->face();
      HalfedgeIter h14 = newHalfedge();
      HalfedgeIter h15 = newHalfedge();
      EdgeIter e6 = newEdge();
      e6->isNew = true;
      FaceIter f2 = newFace();
      h3->next() = h15;
      h3->face() = f2;
      h4->next() = h14;

      h5->face() = f2;
      h10->next() = h4;
      h10->face() = f1;
      h14->twin() = h15;
      h14->next() = h10;
      h14->vertex() = v2;
      h14->edge() = e6;
      h14->face() = f1;
      h15->twin() = h14;
      h15->next() = h5;
      h15->vertex() = v4;
      h15->edge() = e6;
      h15->face() = f2;
      e6->halfedge() = h14;
      f1->halfedge() = h4;
      f2->halfedge() = h3;
    }

    return v4;
  }


  void rearrange(HalfedgeIter h, HalfedgeMesh& mesh, EdgeIter& nextEdge)
  {
    FaceIter f0 = h->face();
    mesh.deleteFace(f0);
    HalfedgeIter next = h->next();
    HalfedgeIter nextNext = next->next();//to be deleted
    HalfedgeIter nextNextTwin = nextNext->twin();//to be deleted
    VertexIter vp = nextNext->vertex();

    HalfedgeIter nextNextTwinNext = nextNextTwin->next();
    vp->halfedge() = nextNextTwinNext;
    next->next() = nextNextTwinNext;
    EdgeIter edge = nextNext->edge();
    if (edge == nextEdge) {
      nextEdge++;
    }
    mesh.deleteEdge(edge);
    mesh.deleteHalfedge(nextNext);
    mesh.deleteHalfedge(nextNextTwin);
  }

  VertexIter HalfedgeMesh::collapseEdge(EdgeIter& e0) {
    // *** Extra Credit ***
    // This method should collapse the given edge and return an iterator to
    // the new vertex created by the collapse.
    EdgeIter nextEdge = e0;
    nextEdge++;

    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->twin();
    VertexIter v0 = h0->vertex();
    VertexIter v1 = h1->vertex();
    int n0 = v0->degree() - 1;
    int n1 = v1->degree() - 1;
    Vector3D midpt = (v1->position + v0->position) / 2;
    bool ish0InTriangle = isInTriangle(h0);
    bool ish1InTriangle = isInTriangle(h1);

    //check connectivity
  //https://stackoverflow.com/questions/27049163/mesh-simplification-edge-collapse-conditions/27049418

    set<VertexCIter> vset0 = v0->getNeighbors();
    set<VertexCIter> vset1 = v1->getNeighbors();
    set<VertexCIter> intersection;

    set_intersection(vset0.begin(), vset0.end(), vset1.begin(), vset1.end(),
      std::inserter(intersection, intersection.begin()));
    int intersectCountShouldBe = (int)ish0InTriangle + (int)ish1InTriangle;
    if (intersection.size() > intersectCountShouldBe) {
      e0 = nextEdge;
      return verticesEnd();
    }

    if (ish0InTriangle) {
      n0--;
    }
    if (ish1InTriangle) {
      n1--;
    }
    int n = n0 + n1;
    //collect e, h, f
    vector<EdgeIter> e(n);
    vector<HalfEdges>h(n);
    vector<FaceIter>f(n);
    HalfedgeIter hs = h1->next();
    for (int i = 0; i < n0; i++) {
      HalfedgeIter ht = hs->twin();
      h[i].h0 = hs;
      h[i].h1 = ht;
      f[i] = ht->face();
      e[i] = hs->edge();
      hs = ht->next();
    }

    hs = h0->next();
    for (int i = n0; i < n; i++) {
      HalfedgeIter ht = hs->twin();
      h[i].h0 = hs;
      h[i].h1 = ht;
      f[i] = ht->face();
      e[i] = hs->edge();
      hs = ht->next();
    }
    //delete vertex, halfedge, face
    deleteVertex(v1);
    if (ish0InTriangle) {
      rearrange(h0, *this, nextEdge);
    }

    if (ish1InTriangle) {
      rearrange(h1, *this, nextEdge);
    }
    //reconnect halfedge->next, halfedge->face,h[i].h0->vertex, face->halfedge
    for (int i = 0; i < n; i++) {
      int inext = next(i, n);
      int iprev = prev(i, n);
      h[i].h1->next() = h[inext].h0;
      h[i].h0->face() = f[iprev];
      h[i].h1->face() = f[i];
      h[i].h0->vertex() = v0;
      f[i]->halfedge() = h[i].h1;
    }

    v0->halfedge() = h[0].h0;
    v0->position = midpt;
    deleteEdge(e0);
    e0 = nextEdge;
    deleteHalfedge(h0);
    deleteHalfedge(h1);
    return v0;
  }

  VertexIter HalfedgeMesh::collapseFace(FaceIter f) {
    // *** Extra Credit ***
    // TODO: (meshEdit)
    // This method should collapse the given face and return an iterator to
    // the new vertex created by the collapse.
    showError("collapseFace() not implemented.");
    return VertexIter();
  }

  FaceIter HalfedgeMesh::eraseVertex(VertexIter v) {
    // *** Extra Credit ***
    // TODO: (meshEdit)
    // This method should replace the given vertex and all its neighboring
    // edges and faces with a single face, returning the new face.

    return FaceIter();
  }

  FaceIter HalfedgeMesh::eraseEdge(EdgeIter e) {
    // *** Extra Credit ***
    // TODO: (meshEdit)
    // This method should erase the given edge and return an iterator to the
    // merged face.

    showError("eraseVertex() not implemented.");
    return FaceIter();
  }

  EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0) {
    // This method should flip the given edge and return an iterator to the
    // flipped edge.
    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h4->next();
    HalfedgeIter h8 = h1->twin();
    HalfedgeIter h9 = h2->twin();
    HalfedgeIter h6 = h4->twin();
    HalfedgeIter h7 = h5->twin();
    // VERTICES
    VertexIter v0 = h2->vertex();
    VertexIter v1 = h0->vertex();
    VertexIter v2 = h6->vertex();
    VertexIter v3 = h3->vertex();
    // EDGES
    EdgeIter e1 = h4->edge();
    EdgeIter e2 = h5->edge();
    EdgeIter e3 = h1->edge();
    EdgeIter e4 = h2->edge();
    // FACES
    FaceIter f0 = h0->face();
    FaceIter f1 = h3->face();

    //reassignment
    h0->next() = h2;
    h0->vertex() = v2;

    h1->next() = h3;
    h1->face() = f1;
    h2->next() = h4;
    h3->next() = h5;
    h3->vertex() = v0;
    h4->next() = h0;
    h4->face() = f0;
    h5->next() = h1;

    v1->halfedge() = h9;
    v3->halfedge() = h1;

    f0->halfedge() = h0;
    f1->halfedge() = h3;

    return e0;
  }

  void HalfedgeMesh::subdivideQuad(bool useCatmullClark) {
    // Unlike the local mesh operations (like bevel or edge flip), we will perform
    // subdivision by splitting *all* faces into quads "simultaneously."  Rather
    // than operating directly on the halfedge data structure (which as you've
    // seen
    // is quite difficult to maintain!) we are going to do something a bit nicer:
    //
    //    1. Create a raw list of vertex positions and faces (rather than a full-
    //       blown halfedge mesh).
    //
    //    2. Build a new halfedge mesh from these lists, replacing the old one.
    //
    // Sometimes rebuilding a data structure from scratch is simpler (and even
    // more
    // efficient) than incrementally modifying the existing one.  These steps are
    // detailed below.

    // Step I: Compute the vertex positions for the subdivided mesh.  Here
    // we're
    // going to do something a little bit strange: since we will have one vertex
    // in
    // the subdivided mesh for each vertex, edge, and face in the original mesh,
    // we
    // can nicely store the new vertex *positions* as attributes on vertices,
    // edges,
    // and faces of the original mesh.  These positions can then be conveniently
    // copied into the new, subdivided mesh.
    if (useCatmullClark) {
      computeCatmullClarkPositions();
    }
    else {
      computeLinearSubdivisionPositions();
    }

    // Step II: Assign a unique index (starting at 0) to each vertex, edge,
    // and
    // face in the original mesh.  These indices will be the indices of the
    // vertices
    // in the new (subdivided mesh).  They do not have to be assigned in any
    // particular
    // order, so long as no index is shared by more than one mesh element, and the
    // total number of indices is equal to V+E+F, i.e., the total number of
    // vertices
    // plus edges plus faces in the original mesh.  Basically we just need a
    // one-to-one
    // mapping between original mesh elements and subdivided mesh vertices.
    assignSubdivisionIndices();

    // Step III: Build a list of quads in the new (subdivided) mesh, as
    // tuples of
    // the element indices defined above.  In other words, each new quad should be
    // of
    // the form (i,j,k,l), where i,j,k and l are four of the indices stored on our
    // original mesh elements.  Note that it is essential to get the orientation
    // right
    // here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces should
    // circulate in the same direction as old faces (think about the right-hand
    // rule).
    vector<vector<Index> > subDFaces;
    vector<Vector3D> subDVertices;
    buildSubdivisionFaceList(subDFaces);
    buildSubdivisionVertexList(subDVertices);

    // Step IV: Pass the list of vertices and quads to a routine that clears
    // the
    // internal data for this halfedge mesh, and builds new halfedge data from
    // scratch,
    // using the two lists.
    rebuild(subDFaces, subDVertices);
  }

  /**
   * Compute new vertex positions for a mesh that splits each polygon
   * into quads (by inserting a vertex at the face midpoint and each
   * of the edge midpoints).  The new vertex positions will be stored
   * in the members Vertex::newPosition, Edge::newPosition, and
   * Face::newPosition.  The values of the positions are based on
   * simple linear interpolation, e.g., the edge midpoints and face
   * centroids.
   */
  void HalfedgeMesh::computeLinearSubdivisionPositions() {
    // For each vertex, assign Vertex::newPosition to
    // its original position, Vertex::position.
    for (VertexIter vertex = verticesBegin(); vertex != verticesEnd(); vertex++) {
      vertex->newPosition = vertex->position;
    }

    // For each edge, assign the midpoint of the two original
    // positions to Edge::newPosition.
    for (EdgeIter edge = edgesBegin(); edge != edgesEnd(); edge++) {
      edge->newPosition = edge->centroid();
    }
    // For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::newPosition.  Note
    // that in general, NOT all faces will be triangles!
    for (FaceIter face = facesBegin(); face != facesEnd(); face++) {
      face->newPosition = face->centroid();
    }
  }

  /**
   * Compute new vertex positions for a mesh that splits each polygon
   * into quads (by inserting a vertex at the face midpoint and each
   * of the edge midpoints).  The new vertex positions will be stored
   * in the members Vertex::newPosition, Edge::newPosition, and
   * Face::newPosition.  The values of the positions are based on
   * the Catmull-Clark rules for subdivision.
   */
  void HalfedgeMesh::computeCatmullClarkPositions() {
    // The implementation for this routine should be
    // a lot like HalfedgeMesh::computeLinearSubdivisionPositions(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules. (These rules are outlined in the Developer Manual.)

    // face
    for (FaceIter face = facesBegin(); face != facesEnd(); face++) {
      face->newPosition = face->centroid();
    }

    // edges
    for (EdgeIter edge = edgesBegin(); edge != edgesEnd(); edge++) {
      HalfedgeIter h = edge->halfedge();
      HalfedgeIter twin = h->twin();
      FaceIter face1 = h->face();
      FaceIter face2 = twin->face();
      edge->newPosition = (face1->newPosition + face2->newPosition) / 2;
    }

    // vertices
    for (VertexIter vertex = verticesBegin(); vertex != verticesEnd(); vertex++) {
      //Q
      HalfedgeIter halfedge = vertex->halfedge();
      HalfedgeIter h = vertex->halfedge();
      int nQ = 0;
      Vector3D Q(0.0, 0.0, 0.0);
      while (true) {
        FaceIter face = h->face();
        Q += face->newPosition;
        nQ++;
        h = h->twin()->next();
        if (h == halfedge) break;
      }
      assert(nQ != 0);
      Q /= nQ;
      //R
      h = vertex->halfedge();
      int nR = 0;
      Vector3D R(0.0, 0.0, 0.0);
      while (true) {
        EdgeIter edge = h->edge();
        R += edge->centroid();
        nR++;
        h = h->twin()->next();
        if (h == halfedge) break;
      }
      assert(nR != 0);
      R /= nR;

      vertex->newPosition = (Q + 2 * R + (nQ - 3) * vertex->position) / nQ;
    }
  }

  /**
   * Assign a unique integer index to each vertex, edge, and face in
   * the mesh, starting at 0 and incrementing by 1 for each element.
   * These indices will be used as the vertex indices for a mesh
   * subdivided using Catmull-Clark (or linear) subdivision.
   */
  void HalfedgeMesh::assignSubdivisionIndices() {
    // Start a counter at zero; if you like, you can use the
    // "Index" type (defined in halfedgeMesh.h)
    Index count = 0;
    // Iterate over vertices, assigning values to Vertex::index
    for (VertexIter vertex = verticesBegin(); vertex != verticesEnd(); vertex++) {
      vertex->index = count++;
    }
    // Iterate over edges, assigning values to Edge::index
    for (EdgeIter edge = edgesBegin(); edge != edgesEnd(); edge++) {
      edge->index = count++;
    }
    // Iterate over faces, assigning values to Face::index
    for (FaceIter face = facesBegin(); face != facesEnd(); face++) {
      face->index = count++;
    }
  }

  /**
   * Build a flat list containing all the vertex positions for a
   * Catmull-Clark (or linear) subdivison of this mesh.  The order of
   * vertex positions in this list must be identical to the order
   * of indices assigned to Vertex::newPosition, Edge::newPosition,
   * and Face::newPosition.
   */
  void HalfedgeMesh::buildSubdivisionVertexList(vector<Vector3D>& subDVertices) {
    // Resize the vertex list so that it can hold all the vertices.
    Index sum = vertices.size() + edges.size() + faces.size();
    subDVertices.reserve(sum);
    // Iterate over vertices, assigning Vertex::newPosition to the
    // appropriate location in the new vertex list.
    for (VertexIter vertex = verticesBegin(); vertex != verticesEnd(); vertex++) {
      subDVertices.push_back(vertex->newPosition);
    }
    // Iterate over edges, assigning Edge::newPosition to the appropriate
    // location in the new vertex list.
    for (EdgeIter edge = edgesBegin(); edge != edgesEnd(); edge++) {
      subDVertices.push_back(edge->newPosition);
    }
    // Iterate over faces, assigning Face::newPosition to the appropriate
    // location in the new vertex list.
    for (FaceIter face = facesBegin(); face != facesEnd(); face++) {
      subDVertices.push_back(face->newPosition);
    }
  }

  /**
   * Build a flat list containing all the quads in a Catmull-Clark
   * (or linear) subdivision of this mesh.  Each quad is specified
   * by a vector of four indices (i,j,k,l), which come from the
   * members Vertex::index, Edge::index, and Face::index.  Note that
   * the ordering of these indices is important because it determines
   * the orientation of the new quads; it is also important to avoid
   * "bowties."  For instance, (l,k,j,i) has the opposite orientation
   * of (i,j,k,l), and if (i,j,k,l) is a proper quad, then (i,k,j,l)
   * will look like a bowtie.
   */
  void HalfedgeMesh::buildSubdivisionFaceList(vector<vector<Index> >& subDFaces) {
    // This routine is perhaps the most tricky step in the construction of
    // a subdivision mesh (second, perhaps, to computing the actual Catmull-Clark
    // vertex positions).  Basically what you want to do is iterate over faces,
    // then for each for each face, append N quads to the list (where N is the
    // degree of the face).  For this routine, it may be more convenient to simply
    // append quads to the end of the list (rather than allocating it ahead of
    // time), though YMMV.  You can of course iterate around a face by starting
    // with its first halfedge and following the "next" pointer until you get
    // back to the beginning.  The tricky part is making sure you grab the right
    // indices in the right order---remember that there are indices on vertices,
    // edges, AND faces of the original mesh.  All of these should get used.  Also
    // remember that you must have FOUR indices per face, since you are making a
    // QUAD mesh!

    // iterate over faces
    for (FaceIter face = facesBegin(); face != facesEnd(); face++) {
      // loop around face
      HalfedgeCIter halfedge = face->halfedge();
      HalfedgeCIter h = face->halfedge();
      Index faceIndex = face->index;
      EdgeCIter edge = h->edge();
      Index edgeIndex = edge->index;
      while (true) {
        // build lists of four indices for each sub-quad
        vector<Index> quad(4);
        quad[0] = faceIndex;
        quad[1] = edgeIndex;
        h = h->next();
        VertexCIter vertex = h->vertex();
        quad[2] = vertex->index;
        edge = h->edge();
        edgeIndex = edge->index;
        quad[3] = edgeIndex;
        // append each list of four indices to face list
        subDFaces.push_back(quad);
        if (h == halfedge) break;
      }
    }
  }



  FaceIter HalfedgeMesh::bevelVertex(VertexIter v) {
    // *** Extra Credit ***
    // This method should replace the vertex v with a face, corresponding to
    // a bevel operation. It should return the new face.  NOTE: This method is
    // responsible for updating the *connectivity* of the mesh only---it does not
    // need to update the vertex positions.  These positions will be updated in
    // HalfedgeMesh::bevelVertexComputeNewPositions (which you also have to
    // implement!)
    int n = v->degree();
    HalfedgeIter hc0_1 = v->halfedge();

    //collect elements
    size_t size = sizeof(HalfEdges) * n;
    HalfEdges* hc = (HalfEdges*)_alloca(size);
    memset(hc, 0, size);
    size = sizeof(EdgeIter) * n;
    EdgeIter* ec = (EdgeIter*)_alloca(size);
    memset(ec, 0, size);
    //    vector<HalfEdges>hc(n);
    size = sizeof(FaceIter) * n;
    FaceIter* fn = (FaceIter*)_alloca(size);
    memset(fn, 0, size);
    //    vector< FaceIter>fn(n);
    HalfedgeIter halfedge = hc0_1;
    for (int i = 0; i < n; i++) {
      hc[i].h1 = halfedge;
      HalfedgeIter halfedgeTwin = halfedge->twin();
      hc[i].h0 = halfedgeTwin;
      ec[i] = hc[i].h0->edge();
      fn[i] = hc[i].h0->face();

      halfedge = halfedgeTwin->next();
    }

    //Allocate new elements
    //halfedges
    size = sizeof(HalfEdges) * n;
    HalfEdges* hi = (HalfEdges*)_alloca(size);
    memset(hi, 0, size);
    //    vector< HalfEdges>hi(n);
    size = sizeof(VertexIter) * n;
    VertexIter* vi = (VertexIter*)_alloca(size);
    memset(vi, 0, size);
    //    vector< VertexIter>vi(n);
    Vector3D pos = v->position;
    size = sizeof(EdgeIter) * n;
    EdgeIter* ei = (EdgeIter*)_alloca(size);
    memset(ei, 0, size);
    //    vector< EdgeIter>ei(n);
    for (int i = 0; i < n; i++) {
      hi[i].h0 = newHalfedge();
      hi[i].h1 = newHalfedge();
      vi[i] = newVertex();
      vi[i]->position = pos;
      ei[i] = newEdge();
    }
    //faces
    FaceIter f = newFace();

    //Reassign
    for (int i = 0; i < n; i++) {
      hc[i].h0->next() = hi[i].h1;
      hc[i].h1->vertex() = vi[i];

      int iNext = next(i, n);
      int iPrev = prev(i, n);
      hi[i].h0->twin() = hi[i].h1;
      hi[i].h0->next() = hi[iPrev].h0;
      hi[i].h0->face() = f;
      hi[i].h0->edge() = ei[i];
      hi[i].h0->vertex() = vi[iNext];

      hi[i].h1->twin() = hi[i].h0;
      hi[i].h1->next() = hc[iNext].h1;
      hi[i].h1->face() = fn[i];
      hi[i].h1->edge() = ei[i];
      hi[i].h1->vertex() = vi[i];

      vi[i]->halfedge() = hi[i].h1;
      ei[i]->halfedge() = hi[i].h1;
    }

    f->halfedge() = hi[0].h0;
    for (int i = 0; i < n; i++) {
      hc[i].h0.HalfedgeIter::~HalfedgeIter();
      hc[i].h1.HalfedgeIter::~HalfedgeIter();
      hi[i].h0.HalfedgeIter::~HalfedgeIter();
      hi[i].h1.HalfedgeIter::~HalfedgeIter();
      ec[i].EdgeIter::~EdgeIter();
      fn[i].FaceIter::~FaceIter();
      vi[i].VertexIter::~VertexIter();
      ei[i].EdgeIter::~EdgeIter();
    }
    return f;
  }

  FaceIter HalfedgeMesh::bevelEdge(EdgeIter e) {
    // *** Extra Credit ***
    // This method should replace the edge e with a face, corresponding to a
    // bevel operation. It should return the new face.  NOTE: This method is
    // responsible for updating the *connectivity* of the mesh only---it does not
    // need to update the vertex positions.  These positions will be updated in
    // HalfedgeMesh::bevelEdgeComputeNewPositions (which you also have to
    // implement!)
    HalfedgeIter h0 = e->halfedge();
    HalfedgeIter h1 = h0->twin();
    VertexIter v0 = h0->vertex();
    VertexIter v1 = h1->vertex();

    int n0 = v0->degree();
    int n1 = v1->degree();
    int n = n0 + n1 - 2;
    //collect em, hm, fm
    vector<HalfEdges> hm(n);
    vector<EdgeIter> em(n);
    vector<FaceIter> fm(n);

    HalfedgeIter hs = h0;
    for (int i = 0; i < n1 - 1; i++) {
      HalfedgeIter next = hs->next();
      HalfedgeIter twin = next->twin();
      hm[i].h0 = next;
      hm[i].h1 = twin;
      em[i] = next->edge();
      fm[i] = hs->face();
      hs = twin;
    }

    hs = h1;
    for (int i = n1 - 1; i < n; i++) {
      HalfedgeIter next = hs->next();
      HalfedgeIter twin = next->twin();
      hm[i].h0 = next;
      hm[i].h1 = twin;
      em[i] = next->edge();
      fm[i] = hs->face();
      hs = twin;
    }

    //allocate new face
    FaceIter f = newFace();

    //edge
    vector<EdgeIter> en(n);
    vector<HalfEdges> hn(n);
    en[0] = e;
    hn[0].h0 = h0;
    hn[0].h1 = h1;
    for (int i = 1; i < n; i++) {
      en[i] = newEdge();
      hn[i].h0 = newHalfedge();
      hn[i].h1 = newHalfedge();
    }

    //vertex
    vector<VertexIter> vn(n);
    vn[0] = v0;
    vn[1] = v1;
    for (int i = 2; i < n; i++) {
      vn[i] = newVertex();
    }

    //reassign
    for (int i = 0; i < n; i++) {
      int iprev = prev(i, n);
      int iprevprev = prev(iprev, n);
      int inext = next(i, n);
      int inextnext = next(inext, n);
      vn[i]->halfedge() = hn[iprev].h0;
      en[i]->halfedge() = hn[i].h0;
      fm[i]->halfedge() = hn[i].h0;

      hm[i].h0->twin() = hm[i].h1;
      hm[i].h1->twin() = hm[i].h0;
      hn[i].h0->twin() = hn[i].h1;
      hn[i].h1->twin() = hn[i].h0;

      hm[i].h1->next() = hn[inext].h0;
      hn[i].h0->next() = hm[i].h0;
      hn[i].h1->next() = hn[iprev].h1;

      hm[i].h0->vertex() = vn[inextnext];
      hn[i].h0->vertex() = vn[inext];
      hn[i].h1->vertex() = vn[inextnext];

      hn[i].h0->edge() = en[i];
      hn[i].h1->edge() = en[i];

      hn[i].h0->face() = fm[i];
      hn[i].h1->face() = f;
    }

    for (int i = 0; i < n1 - 1; i++) {
      hm[i].h0->vertex()->position = v1->position;
    }
    for (int i = n1 - 1; i < n; i++) {
      hm[i].h0->vertex()->position = v0->position;
    }

    f->halfedge() = hn[0].h1;
    return f;
  }

  FaceIter HalfedgeMesh::bevelFace(FaceIter f) {
    // *** Extra Credit ***
    // This method should replace the face f with an additional, inset face
    // (and ring of faces around it), corresponding to a bevel operation. It
    // should return the new face.  NOTE: This method is responsible for updating
    // the *connectivity* of the mesh only---it does not need to update the vertex
    // positions.  These positions will be updated in
    // HalfedgeMesh::bevelFaceComputeNewPositions (which you also have to
    // implement!)
    int n = f->degree();
    HalfedgeIter h0_0 = f->halfedge();

    //collect elements
    size_t size = sizeof(VertexIter) * n;
    VertexIter* v = (VertexIter*)_alloca(size);
    memset(v, 0, size);
    HalfedgeIter halfedge = h0_0;
    for (int i = 0; i < n; i++) {
      VertexIter vertex = halfedge->vertex();
      v[i] = vertex;
      halfedge = halfedge->next();
    }

    size = sizeof(EdgeIter) * n;
    EdgeIter* e = (EdgeIter*)_alloca(size);
    memset(e, 0, size);
    halfedge = h0_0;
    for (int i = 0; i < n; i++) {
      e[i] = halfedge->edge();
      halfedge = halfedge->next();
    }

    size = sizeof(HalfEdges) * n;
    HalfEdges* h = (HalfEdges*)_alloca(size);
    memset(h, 0, size);
    halfedge = h0_0;
    for (int i = 0; i < n; i++) {
      HalfEdges& hi = h[i];
      hi.h0 = halfedge;
      hi.h1 = halfedge->twin();
      halfedge = halfedge->next();
    }

    //allocate new elements
    size = sizeof(VertexIter) * n;
    VertexIter* vi = (VertexIter*)_alloca(size);
    memset(vi, 0, size);
    for (int i = 0; i < n; i++) {
      vi[i] = newVertex();
      vi[i]->position = v[i]->position;
    }

    size = sizeof(EdgeIter) * n;
    EdgeIter* ei = (EdgeIter*)_alloca(size);
    memset(ei, 0, size);
    for (int i = 0; i < n; i++) {
      ei[i] = newEdge();
    }

    size = sizeof(EdgeIter) * n;
    EdgeIter* ec = (EdgeIter*)_alloca(size);
    memset(ec, 0, size);
    for (int i = 0; i < n; i++) {
      ec[i] = newEdge();
    }

    size = sizeof(FaceIter) * n;
    FaceIter* fn = (FaceIter*)_alloca(size);
    memset(fn, 0, size);
    for (int i = 0; i < n; i++) {
      fn[i] = newFace();
    }

    size = sizeof(HalfEdges) * n;
    HalfEdges* hi = (HalfEdges*)_alloca(size);
    memset(hi, 0, size);
    for (int i = 0; i < n; i++) {
      hi[i].h0 = newHalfedge();
      hi[i].h1 = newHalfedge();
    }

    size = sizeof(HalfEdges) * n;
    HalfEdges* hc = (HalfEdges*)_alloca(size);
    memset(hc, 0, size);
    for (int i = 0; i < n; i++) {
      hc[i].h0 = newHalfedge();
      hc[i].h1 = newHalfedge();
    }
    //reassign elements
    for (int i = 0; i < n; i++) {
      HalfedgeIter halfedge = h[i].h0;
      halfedge->next() = hc[next(i, n)].h0;
      halfedge->face() = fn[i];
    }
    //hc1_0
    for (int i = 0; i < n; i++) {
      HalfedgeIter halfedge = hc[i].h0;
      int iPrev = prev(i, n);
      halfedge->next() = hi[iPrev].h1;
      halfedge->twin() = hc[i].h1;
      halfedge->vertex() = v[i];
      halfedge->edge() = ec[i];
      halfedge->face() = fn[iPrev];
    }

    //hc1_1
    for (int i = 0; i < n; i++) {
      HalfedgeIter halfedge = hc[i].h1;
      halfedge->next() = h[i].h0;
      halfedge->twin() = hc[i].h0;
      halfedge->vertex() = vi[i];
      halfedge->edge() = ec[i];
      halfedge->face() = fn[i];
    }

    //hi0_1
    for (int i = 0; i < n; i++) {
      HalfedgeIter halfedge = hi[i].h1;

      halfedge->next() = hc[i].h1;
      halfedge->twin() = hi[i].h0;
      int iNext = next(i, n);
      halfedge->vertex() = vi[iNext];
      halfedge->edge() = ei[i];
      halfedge->face() = fn[i];
    }

    //hi0_0
    for (int i = 0; i < n; i++) {
      HalfedgeIter halfedge = hi[i].h0;
      int iNext = next(i, n);
      halfedge->next() = hi[iNext].h0;
      halfedge->twin() = hi[i].h1;

      halfedge->vertex() = vi[i];
      halfedge->edge() = ei[i];
      halfedge->face() = f;
    }
    // VERTICES
    for (int i = 0; i < n; i++) {
      VertexIter vertex = vi[i];
      vertex->halfedge() = hi[i].h0;
    }

    // EDGES
    for (int i = 0; i < n; i++) {
      EdgeIter edge = ec[i];
      edge->halfedge() = hc[i].h0;
    }
    for (int i = 0; i < n; i++) {
      EdgeIter edge = ei[i];
      edge->halfedge() = hi[i].h0;
    }
    // FACES
    for (int i = 0; i < n; i++) {
      FaceIter face = fn[i];
      face->halfedge() = h[i].h0;
    }
    f->halfedge() = hi[0].h0;

    for (int i = 0; i < n; i++) {
      h[i].h0.HalfedgeIter::~HalfedgeIter();
      h[i].h1.HalfedgeIter::~HalfedgeIter();
      hi[i].h0.HalfedgeIter::~HalfedgeIter();
      hi[i].h1.HalfedgeIter::~HalfedgeIter();
      hc[i].h0.HalfedgeIter::~HalfedgeIter();
      hc[i].h1.HalfedgeIter::~HalfedgeIter();
      ei[i].EdgeIter::~EdgeIter();
      ec[i].EdgeIter::~EdgeIter();
      fn[i].FaceIter::~FaceIter();
      v[i].VertexIter::~VertexIter();
      vi[i].VertexIter::~VertexIter();
      e[i].EdgeIter::~EdgeIter();
    }
    return f;
  }

  void HalfedgeMesh::bevelFaceComputeNewPositions(
    vector<Vector3D>& originalVertexPositions,
    vector<HalfedgeIter>& newHalfedges, double normalShift,
    double tangentialInset) {
    // *** Extra Credit ***
    // Compute new vertex positions for the vertices of the beveled face.
    //
    // These vertices can be accessed via newHalfedges[i]->vertex()->position for
    // i = 1, ..., newHalfedges.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the originalVertexPositions array) to compute an offset vertex
    // position.
    //
    // Note that there is a 1-to-1 correspondence between halfedges in
    // newHalfedges and vertex positions
    // in orig.  So, you can write loops of the form
    //
    // for( int i = 0; i < newHalfedges.size(); i++ )
    // {
    //    Vector3D pi = originalVertexPositions[i]; // get the original vertex
    //    position correponding to vertex i
    // }
    //
    int n = newHalfedges.size();
    vector<Vector3D> normals(n);
    //in case originalVertexPositions's 3D points are not in a plane, I use bisector to calculate new shift points
    for (int i = 0; i < n; i++)
    {
      Vector3D pi = originalVertexPositions[i];
      int inext = next(i, n);
      Vector3D pNext = originalVertexPositions[inext];
      normals[i] = (pNext - pi);
      normals[i].normalize();
    }

    for (int i = 0; i < n; i++)
    {
      int iprev = prev(i, n);
      const Vector3D& nomalPrev = normals.at(iprev);
      const Vector3D& nomal = normals.at(i);
      Vector3D bisector = nomal - nomalPrev;
      bisector.normalize();
      double cosa = dot(bisector, nomal);
      double sina = sqrt(1 - cosa * cosa);
      Vector3D up = cross(nomalPrev, nomal);
      up.normalize();
      newHalfedges[i]->vertex()->position += (bisector * normalShift * sina + up * tangentialInset);
    }
  }

  void HalfedgeMesh::bevelVertexComputeNewPositions(
    Vector3D originalVertexPosition, vector<HalfedgeIter>& newHalfedges,
    double tangentialInset) {
    // *** Extra Credit ***
    // Compute new vertex positions for the vertices of the beveled vertex.
    //
    // These vertices can be accessed via newHalfedges[i]->vertex()->position for
    // i = 1, ..., hs.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the orig array) to compute an offset vertex position.
    for (auto& halfedge : newHalfedges) {
      Vector3D pos = halfedge->twin()->vertex()->position;

      //Vector3D pos = halfedge->vertex()->position;
      Vector3D vec = originalVertexPosition - pos;
      double length = vec.norm() * tangentialInset;
      Vector3D dir = vec.unit() * length;
      halfedge->vertex()->position += dir;
    }

  }

  void HalfedgeMesh::bevelEdgeComputeNewPositions(
    vector<Vector3D>& originalVertexPositions,
    vector<HalfedgeIter>& newHalfedges, double tangentialInset) {
    // *** Extra Credit ***
    // Compute new vertex positions for the vertices of the beveled edge.
    //
    // These vertices can be accessed via newHalfedges[i]->vertex()->position for
    // i = 1, ..., newHalfedges.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the orig array) to compute an offset vertex position.
    //
    // Note that there is a 1-to-1 correspondence between halfedges in
    // newHalfedges and vertex positions
    // in orig.  So, you can write loops of the form

    for (int i = 0; i < newHalfedges.size(); i++) {
      HalfedgeIter halfedge = newHalfedges[i];
      Vector3D pos = halfedge->twin()->vertex()->position;
      Vector3D pi = originalVertexPositions[i];
      //Vector3D pos = halfedge->vertex()->position;
      Vector3D vec = pi - pos;
      double length = vec.norm() * tangentialInset;
      Vector3D dir = vec.unit() * length;
      halfedge->vertex()->position += dir;
    }

  }

  void HalfedgeMesh::splitPolygons(vector<FaceIter>& fcs) {
    for (auto f : fcs) splitPolygon(f);
  }

  void HalfedgeMesh::splitPolygon(FaceIter f) {
    // *** Extra Credit ***
    // TODO: (meshedit) 
    // Triangulate a polygonal face
    showError("splitPolygon() not implemented.");
  }

  EdgeRecord::EdgeRecord(EdgeIter& _edge) : edge(_edge) {
    // *** Extra Credit ***
    // (meshEdit)
    // Compute the combined quadric from the edge endpoints.
    // -> Build the 3x3 linear system whose solution minimizes the quadric error
    //    associated with these two endpoints.
    // -> Use this system to solve for the optimal position, and store it in
    //    EdgeRecord::optimalPoint.
    // -> Also store the cost associated with collapsing this edg in
    //    EdgeRecord::Cost.
    Matrix4x4 quadric = edge->getQuadric();
    Matrix3x3 A;
    A.column(0) = quadric.column(0).to3D();
    A.column(1) = quadric.column(1).to3D();
    A.column(2) = quadric.column(2).to3D();
    Vector3D b = -quadric.column(3).to3D();
    optimalPoint = A.inv() * b;
    Vector4D optimalPoint4D(optimalPoint, 1);
    this->score = dot(optimalPoint4D, (quadric * optimalPoint4D));
  }

  void MeshResampler::upsample(HalfedgeMesh& mesh)
    // This routine should increase the number of triangles in the mesh using Loop
    // subdivision.
  {
    // Compute new positions for all the vertices in the input mesh, using
    // the Loop subdivision rule, and store them in Vertex::newPosition.
    // -> At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh.
    for (VertexIter vertex = mesh.verticesBegin(); vertex != mesh.verticesEnd(); vertex++) {
      if (!vertex->isBoundary()) {
        std::pair<CS248::Vector3D, int> sumD = vertex->neighborhoodSumAndDegree();
        Vector3D sum = sumD.first;
        int n = sumD.second;
        double u = 0;
        if (n == 3) {
          u = 3.0 / 16;
        }
        else {
          u = 3.0 / (8 * n);
        }
        vertex->newPosition = (1 - n * u) * vertex->position + u * sum;
      }
      vertex->isNew = false;
    }
    for (EdgeIter edge = mesh.edgesBegin(); edge != mesh.edgesEnd(); edge++) {
      edge->isNew = false;
    }
    //For boundary edges we need special treatment
    for (FaceIter boundary = mesh.boundariesBegin(); boundary != mesh.boundariesEnd(); boundary++) {
      HalfedgeIter h0 = boundary->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter _h2 = h1->next();
      HalfedgeIter h2 = _h2;
      VertexIter v0 = h0->vertex();
      VertexIter v1 = h1->vertex();
      VertexIter v2 = h2->vertex();
      while (true) {
        v1->newPosition = 1.0 / 8 * v0->position + 3.0 / 4 * v1->position + 1.0 / 8 * v2->position;
        h2 = h2->next();
        if (h2 == _h2) break;
        v0 = v1;
        v1 = v2;
        v2 = h2->vertex();
      }
    }
    // -> Next, compute the updated vertex positions associated with edges, and
    //    store it in Edge::newPosition.
    for (EdgeIter edge = mesh.edgesBegin(); edge != mesh.edgesEnd(); edge++) {
      Vector3D center = edge->centroid();
      if (edge->isBoundary()) {
        edge->newPosition = center;
      }
      else {
        HalfedgeCIter h0 = edge->halfedge();
        HalfedgeCIter h1 = h0->twin();
        HalfedgeCIter h0_back = h0->next()->next();
        HalfedgeCIter h1_back = h1->next()->next();
        Vector3D p0 = h0_back->vertex()->position;
        Vector3D p1 = h1_back->vertex()->position;
        edge->newPosition = 3.0 / 4 * center + 1.0 / 8 * (p0 + p1);
      }

    }

    // -> Next, we're going to split every edge in the mesh, in any order.  For
    //    future reference, we're also going to store some information about which
    //    subdivided edges come from splitting an edge in the original mesh, and
    //    which edges are new, by setting the flat Edge::isNew. Note that in this
    //    loop, we only want to iterate over edges of the original mesh.
    //    Otherwise, we'll end up splitting edges that we just split (and the
    //    loop will never end!)
    EdgeIter edge = mesh.edgesBegin();
    EdgeIter end = mesh.edgesEnd();
    EdgeIter lastEdge = --end;
    while (true) {
      VertexIter vertex = mesh.splitEdge(edge);
      vertex->newPosition = edge->newPosition;
      if (edge == lastEdge) break;
      edge++;
    }

    // -> Now flip any new edge that connects an old and new vertex.
    for (EdgeIter edge = mesh.edgesBegin(); edge != mesh.edgesEnd(); edge++) {
      if (edge->isNew) {
        HalfedgeCIter h0 = edge->halfedge();
        HalfedgeCIter h1 = h0->twin();
        if (!h0->vertex()->isNew || !h1->vertex()->isNew) {
          mesh.flipEdge(edge);
        }
      }
    }

    // -> Finally, copy the new vertex positions into final Vertex::position.
    for (VertexIter vertex = mesh.verticesBegin(); vertex != mesh.verticesEnd(); vertex++) {
      vertex->position = vertex->newPosition;
    }
    // Each vertex and edge of the original surface can be associated with a
    // vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to
    // *first* compute the new positions
    // using the connectity of the original (coarse) mesh; navigating this mesh
    // will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.  We
    // will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.
  }

  static void removeVertexEdges(VertexIter v, set<EdgeRecord>& toBeRemoved)
  {
    HalfedgeIter _h = v->halfedge();
    HalfedgeIter h = _h;
    while (true) {
      EdgeIter edge = h->edge();
      toBeRemoved.insert(edge->record);
      h = h->twin()->next();
      if (h == _h)break;
    }
  }

  void MeshResampler::downsample(HalfedgeMesh& mesh) {
    // *** Extra Credit ***
    // (meshEdit)
    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in Face::quadric
    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in Vertex::quadric
    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an EdgeRecord for each edge and sticking it in the
    //    queue.
    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.
    for (FaceIter it = mesh.facesBegin(); it != mesh.facesEnd(); it++) {
      it->calcQuadric();
    }
    for (VertexIter it = mesh.verticesBegin(); it != mesh.verticesEnd(); it++) {
      it->calcQuadric();
    }
    MutablePriorityQueue<EdgeRecord> queue;
    for (EdgeIter it = mesh.edgesBegin(); it != mesh.edgesEnd(); it++) {
      EdgeRecord record(it);
      it->record = record;
      queue.insert(record);
    }
    int nFaces = mesh.nFaces();
    int targetFaces = nFaces / 4;
    while (true) {
      int n = mesh.nFaces();
      if (n <= targetFaces) break;
      EdgeRecord record = queue.top();
      queue.pop();
      EdgeIter edge = record.edge;

      HalfedgeIter h0 = edge->halfedge();
      HalfedgeIter h1 = h0->twin();
      VertexIter v0 = h0->vertex();
      VertexIter v1 = h1->vertex();
      Matrix4x4 quatric = v0->quadric + v1->quadric;
      set<EdgeRecord> toBeRemoved;
      removeVertexEdges(v0, toBeRemoved);
      removeVertexEdges(v1, toBeRemoved);
      toBeRemoved.erase(record);
      for (auto& record : toBeRemoved) {
        if (!queue.remove(record)) {
          toBeRemoved.erase(record);
        }
      }
      VertexIter vertex = mesh.collapseEdge(edge);
      if (vertex != mesh.verticesEnd()) {
        vertex->position = record.optimalPoint;
        vertex->quadric = quatric;
        HalfedgeIter _h = vertex->halfedge();
        HalfedgeIter h = _h;
        while (true) {

          EdgeIter edge = h->edge();
          EdgeRecord record(edge);
          edge->record = record;
          queue.insert(record);
          h = h->twin()->next();
          if (h == _h) break;
        }
      }
      else {
        for (auto& record : toBeRemoved) {
          queue.insert(record);
        }

      }
    }



  }

  void MeshResampler::resample(HalfedgeMesh& mesh) {
    // *** Extra Credit ***
    // (meshEdit)
    // Compute the mean edge length.
    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions
    double sum = 0.0;
    for (auto edge = mesh.edgesBegin(); edge != mesh.edgesEnd(); edge++) {
      sum += edge->length();
    }
    double meanLength = sum / mesh.nEdges();
    double upperLimit = meanLength * 4 / 3;
    double lowerLimit = meanLength * 4 / 5;

    for (int i = 0; i < 5; i++) {
      for (auto edge = mesh.edgesBegin(); edge != mesh.edgesEnd(); edge++) {
        double length = edge->length();
        if (length > upperLimit) {
          mesh.splitEdge(edge);
        }
      }


      for (auto edge = mesh.edgesBegin(); edge != mesh.edgesEnd();) {
        double length = edge->length();
        if (length < lowerLimit) {
          mesh.collapseEdge(edge);
        }
        else {
          edge++;
        }
      }
      //flip
      for (auto edge = mesh.edgesBegin(); edge != mesh.edgesEnd();) {
        auto oldEdge = edge++;
        auto h0 = oldEdge->halfedge();
        auto h1 = h0->twin();
        auto v0 = h0->vertex();
        auto v1 = h1->vertex();
        int a0 = v0->degree();
        int a1 = v1->degree();
        auto ch0 = h0->next()->next();
        auto ch1 = h1->next()->next();
        auto cv0 = ch0->vertex();
        auto cv1 = ch1->vertex();
        int b0 = cv0->degree();
        int b1 = cv1->degree();
        int div0 = abs(a0 - 6) + abs(a1 - 6) + abs(b0 - 6) + abs(b1 - 6);
        int div1 = abs(a0 - 1 - 6) + abs(a1 - 1 - 6) + abs(b0 + 1 - 6) + abs(b1 + 1 - 6);
        if (div1 < div0) {
          mesh.flipEdge(oldEdge);
        }
      }
      //computer new vertex position
      for (auto vertex = mesh.verticesBegin(); vertex != mesh.verticesEnd(); vertex++) {
        Vector3D c = vertex->neighborhoodCentroid();
        Vector3D p = vertex->position;
        Vector3D v = c - p;
        Vector3D N = vertex->normal();
        v = v - dot(N, v) * N;
        vertex->newPosition = p + 1.0 / 5 * v;
      }
      for (auto vertex = mesh.verticesBegin(); vertex != mesh.verticesEnd(); vertex++) {
        vertex->position = vertex->newPosition;
      }

    }
  }



}  // namespace CS248
