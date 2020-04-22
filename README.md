<p><a href="#2a84ce044b53bb3902f080782d296b09">Vector2() - Vectors in 2 dimensional euclidean space.</a>

<h1 id="2a84ce044b53bb3902f080782d296b09">Class: Vector2() <h1>
<p>Vectors in 2 dimensional euclidean space.</p>

<h2>Normalize</h2>
Normalize a copy of vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to clone and normalize

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).Normalize() == x * 0.6 + y * 0.8)
</pre>
<hr>
<h2>R180</h2>
Rotate  a copy of a vector by 180 degrees.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to clone and rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.R180() == -x)
</pre>
<hr>
<h2>R270</h2>
Rotate a vector by 270 degrees anticlockwise.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to clone and rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.R270() == -y)
</pre>
<hr>
<h2>R90</h2>
Rotate a copy of a vector by 90 degrees anticlockwise.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to clone and rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.R90() == y)
</pre>
<hr>
<h2>Swap</h2>
Swap the components of a copy of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to clone and swap

</table>

<h4>Examples</h4>
<pre language="python">
z = x + y * 2
Z = z.Swap()
assert(z != Z)
</pre>
<hr>
<h2>__abs__</h2>
Length of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector whose length is required

</table>

<h4>Examples</h4>
<pre language="python">
assert(abs(x * 3 + y * 4) == 5)
</pre>
<hr>
<h2>__add__ **+ **</h2>
Add the second vector to a copy of the first vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned and added to

<tr><td><b><i>other</i></b><td>Vector to add

</table>

<h4>Examples</h4>
<pre language="python">
assert(x + y == Vector2(1, 1))
</pre>
<hr>
<h2>__iadd__ **+=**</h2>
Add the second vector to the first vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be added to

<tr><td><b><i>other</i></b><td>Vector to add

</table>

<h4>Examples</h4>
<pre language="python">
X = x
X = X + x
assert(x == Vector2(1, 0) and X == Vector2(2, 0))
</pre>
<hr>
<h2>__imul__ ***=**</h2>
Multiply a vector by a scalar.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be multiplied

<tr><td><b><i>scalar</i></b><td>Scalar to multiply by

</table>

<h4>Examples</h4>
<pre language="python">
X = x
X = X * 2
assert(x == Vector2(1, 0) and X == Vector2(2, 0))
</pre>
<hr>
<h2>__init__</h2>
Create a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>New vector to be initialized

<tr><td><b><i>x = 0</i></b><td>X coordinate

<tr><td><b><i>y = 0</i></b><td>Y coordinate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x + y == Vector2(1, 1))
</pre>
<hr>
<h2>__isub__ **-=**</h2>
Subtract the second vector from the first vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be subtracted from

<tr><td><b><i>other</i></b><td>Vector to subtract

</table>

<h4>Examples</h4>
<pre language="python">
X = x
X = X - x
assert(x == Vector2(1, 0) and X == Vector2(0, 0))
</pre>
<hr>
<h2>__itruediv__ **/=**</h2>
Divide a vector by a scalar.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be divided

<tr><td><b><i>scalar</i></b><td>Scalar to divide by

</table>

<h4>Examples</h4>
<pre language="python">
X = x * 4
X = X / 2
assert(x == Vector2(1, 0) and X == Vector2(2, 0))
</pre>
<hr>
<h2>__mul__ *** **</h2>
Multiply a copy of vector by a scalar.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned and multiplied

<tr><td><b><i>scalar</i></b><td>Scalar to multiply by

</table>

<h4>Examples</h4>
<pre language="python">
assert(x * 2 + y * 3 == Vector2(2, 3))
</pre>
<hr>
<h2>__neg__ **- **</h2>
Rotate a copy of a vector by 180 degrees.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned and negated

</table>

<h4>Examples</h4>
<pre language="python">
assert(- y == Vector2(0, -1))
</pre>
<hr>
<h2>__repr__</h2>
String representation of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be represented

</table>

<h4>Examples</h4>
<pre language="python">
assert(repr(x + y * 2) == 'Vector2(1, 2)')
</pre>
<hr>
<h2>__sub__ **- **</h2>
Subtract the second vector from a copy of the first vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned and subtracted from

<tr><td><b><i>other</i></b><td>Vector to subtract

</table>

<h4>Examples</h4>
<pre language="python">
assert(x - y == Vector2(1, -1))
</pre>
<hr>
<h2>__truediv__ **/ **</h2>
Divide a copy of a vector by a scalar.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned and divided

<tr><td><b><i>scalar</i></b><td>Scalar to divide by

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 4 + y * 2) / 2 == x * 2 + y)
</pre>
<hr>
<h2>angle</h2>
Angle in radians anticlockwise that the first vector must be rotated to point along the second vector normalized to the range: -pi to +pi.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>o</i></b><td>Left hand vector we are looking along

<tr><td><b><i>p</i></b><td>Right hand vector we want anticlockwise angle to

</table>

<h4>Examples</h4>
<pre language="python">
for i in range(-179, +179): # Anticlockwise angle from x
assert(Vector2.close
(x.angle(x * math.cos(dr(i)) + y * math.sin(dr(i))), dr(i)))
</pre>
<hr>
<h2>area</h2>
Signed area of the parallelogram defined by the two vectors. The area is negative if the second vector appears to the right of the first if they are both placed at the origin and the observer stands against the z-axis in a left handed coordinate system.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Left hand vector

<tr><td><b><i>other</i></b><td>Right hand vector

</table>

<h4>Examples</h4>
<pre language="python">
assert((x + y).area(-x + y) == 2)
</pre>
<hr>
<h2>clone</h2>
Clone a vector to allow it to be modified by other operations.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned

</table>

<h4>Examples</h4>
<pre language="python">
z = x + y * 2
Z = z.clone()
assert(z == Z)
</pre>
<hr>
<h2>close</h2>
Whether two numbers are close.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>a</i></b><td>First number to compare

<tr><td><b><i>b</i></b><td>Second number to compare

</table>

<h4>Examples</h4>
<pre language="python">
assert(Vector2.close(0, Vector2.closeRadius() / 2))
</pre>
<hr>
<h2>cos</h2>
cos(angle between two vectors).
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>First vector

<tr><td><b><i>other</i></b><td>Second vector

</table>

<h4>Examples</h4>
<pre language="python">
assert(Vector2.close((x + y).cos(y), 1 / r2))
assert(Vector2.close( x.cos(x + yr3), 0.5))
</pre>
<hr>
<h2>distance</h2>
Distance between the points identified by two vectors when placed on the same point.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to start point

<tr><td><b><i>other</i></b><td>Vector to end point

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).distance (-(x * 3 + y * 4)) == 10)
</pre>
<hr>
<h2>distance2</h2>
Distance squared between the points identified

       by two vectors when placed on the same point.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to start point

<tr><td><b><i>other</i></b><td>Vector to end point

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).distance2(-(x * 3 + y * 4)) == 100)
</pre>
<hr>
<h2>dot</h2>
Dot product of two vectors.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>First vector

<tr><td><b><i>other</i></b><td>Second vector

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 2 + y).dot(x + y * 3) == 5)
</pre>
<hr>
<h2>length</h2>
Length of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector whose length is required

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).length() == 5)
</pre>
<hr>
<h2>length2</h2>
Length squared of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vectors whose length squared is required

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).length2() == 25)
</pre>
<hr>
<h2>normalize</h2>
Normalize a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to normalize

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).clone().normalize() == x * 0.6 + y * 0.8)
</pre>
<hr>
<h2>r180</h2>
Rotate a vector by 180 degrees.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.clone().r180() == -x)
</pre>
<hr>
<h2>r270</h2>
Rotate a copy of a vector by 270 degrees anticlockwise.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.clone().r270() == -y)
</pre>
<hr>
<h2>r90</h2>
Rotate a vector by 90 degrees anticlockwise.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.clone().r90() == y)
</pre>
<hr>
<h2>sin</h2>
sin(angle between two vectors).
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Left hand vector

<tr><td><b><i>other</i></b><td>Right hand vector

</table>

<h4>Examples</h4>
<pre language="python">
assert(Vector2.close((x + y).sin(y), 1 / r2))
assert(Vector2.close( x.sin(x + yr3), r3 / 2))
</pre>
<hr>
<h2>smallestAngleToNormalPlane</h2>
The smallest angle between the second vector and a plane normal to the first vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>a</i></b><td>Vector normal to plane

<tr><td><b><i>b</i></b><td>Vector at angle to plane

</table>

<h4>Examples</h4>
<pre language="python">
assert(Vector2.close(dr( 0), y.smallestAngleToNormalPlane( x))) # First vector is y, second vector is 0 degrees anti-clockwise from x axis
assert(Vector2.close(dr(+45), y.smallestAngleToNormalPlane( x + y))) # +45
assert(Vector2.close(dr(+90), y.smallestAngleToNormalPlane( y))) # +90
assert(Vector2.close(dr(+45), y.smallestAngleToNormalPlane(-x + -y))) # +135
assert(Vector2.close(dr( 0), y.smallestAngleToNormalPlane(-x))) # +180
assert(Vector2.close(dr(+45), y.smallestAngleToNormalPlane(-x + -y))) # +225
assert(Vector2.close(dr(+90), y.smallestAngleToNormalPlane( -y))) # +270
assert(Vector2.close(dr(+45), y.smallestAngleToNormalPlane(-x + -y))) # +315
assert(Vector2.close(dr( 0), y.smallestAngleToNormalPlane( x))) # +360
</pre>
<hr>
<h2>swap</h2>
Swap the components of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to swap

</table>

<h4>Examples</h4>
<pre language="python">
assert((x + y * 2).swap() == y + x * 2)
</pre>
<hr>
