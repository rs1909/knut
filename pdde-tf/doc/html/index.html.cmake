<?xml version="1.0" encoding="iso-8859-1"?>
<!DOCTYPE html
    PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
	<title>@PACKAGE_NAME@</title>
	<style>h3.fn,span.fn { margin-left: 1cm; text-indent: -1cm; }
		a:link { color: #2B5E82; text-decoration: none }
		a:visited { color: #672967; text-decoration: none }
		td.postheader { font-family: sans-serif }
		tr.address { font-family: sans-serif }
		body { background: "#ffffff"; color: black; }
	</style>
</head>

<body>
<table width="640" align="center"> <tr><td>
<p>
<table cellpadding="2" cellspacing="1" border="0" width="100%" bgcolor="#e5e5e5">
	<tr>
		<th bgcolor="#333366" style="color: white"> <h2>@PACKAGE_NAME@: A continuation and bifurcation software for delay-differential equations</h2></th>
	</tr>
</table>

<p />
<b>Note</b>
<hr size="1" width="100%" color="black" />
<p>
This software is a new and renamed version of <a href="http://seis.bris.ac.uk/~rs1909/pdde/">PDDE-CONT</a>.
</p>

<b>Capabilities</b>
<hr size="1" width="100%" color="black" />
<p>
@PACKAGE_NAME@ can continue periodic solutions and periodic solution bifurcations of autonomous and periodically forced delay-differential equations that are in the general form</p>
<p align="center">
<i>dx(t)/dt = f(t, x(t-&tau;<sub>0</sub>), x(t-&tau;<sub>1</sub>), ..., x(t-&tau;<sub>r</sub>), &lambda;).</i>
</p>
<p>
To start the continuation, the software needs an initial periodic solution <i>x<sub>0</sub>(t)</i> at some parameter value <i>&lambda;<sub>0</sub></i>. This starting periodic
solution has to be specified either in the system definition (C++ source file) or it can be loaded from
an input file. The periodic solutions can bifurcate in several ways. If one of the three common codimension-one
bifurcations (fold, period doubling, Neimark-Sacker) is found along the branch of periodic
solutions, the point can be used as a starting point for continuing the branch of bifurcation points in a two parameter space.
</p>

<p>
Equilibria of autonomous systems can be also handled, as periodic solutions of periodic
systems by keeping the period constant. Although it is very inefficient, this capability might be useful,
when looking for periodic orbits arising at Hopf bifurcations points. Hopf bifurcations of fixed points
are detected as Neimark-Sacker bifurcations, and the periodic solutions emanating from such points can
be continued by switching to the periodic solution branch.
</p>

<b> Documentation </b>
<hr size="1" width="100%" color="black" />
<p>
See the <a href="manual.pdf">User's Manual</a>.
</p>

<b> Download </b>
<hr size="1" width="100%" color="black" />

<p>
<table cellpadding="2" cellspacing="1" border="0" width="100%" bgcolor="#F6F6F6">
<tr style="color: white"><th bgcolor="#333366">Current version (@PACKAGE_VERSION@)</th> <th bgcolor="#333366"> Package </th>
</tr>
<tr>
<td> Source code </td> <td> <a href="@PACKAGE_NAME@-@PACKAGE_VERSION@.tar.bz2"> @PACKAGE_NAME@-@PACKAGE_VERSION@.tar.bz2 </a></td>
</tr>
<tr>
<td> Mac Application (Intel x86_64) </td> <td> <a href="@PACKAGE_NAME@-@PACKAGE_VERSION@.dmg"> @PACKAGE_NAME@-@PACKAGE_VERSION@.dmg</a></td></tr>
</table>
</p><p/>

<b>Installation from the source code</b>
<hr size="1" width="100%" color="black" /><p>
The software can be built using the cross-platform build facility <a href="http://www.cmake.org">CMake</a>. The build also requires a C and C++ compiler, preferably <a href="http://gcc.gnu.org">GCC</a>, and the <a href="http://math-atlas.sourceforge.net/">ATLAS</a> linear algebra library. The (optional) graphical user interface requires the <a href="http://www.qtsoftware.com/products/qt/downloads">Qt</a> GUI toolkit.
</p>

<b>License</b>
<hr size="1" width="100%" color="black" />
<p>
	The software is available free of charge. Copying of the software is governed by the GNU General Public License <a href="http://www.gnu.org/licenses/gpl.html">(GPL)</a> version 2.
	In addition to the GPL, the author of any scientific or non-scientific publication who used the software is required to refer to the User's Manual in the publication where the software or its modified version was used <a href="#man">[1]</a>.
</p>

<b>
Third party Licenses used in KNUT</b>
<hr size="1" width="100%" color="black" /><p>
<a href="http://www.cise.ufl.edu/research/sparse/umfpack">UMFPACK</a> and <a href="http://www.cise.ufl.edu/research/sparse/amd">AMD</a> 
can be distributed under the GNU Lesser General Public License <a href="http://www.gnu.org/licenses/lgpl.html">(LGPL)</a>
</p>

<b> Related publications </b>
<hr size="1" width="100%" color="black" />
<p>
<a name="man">[1]</a> R. Szalai, <i> Knut: A continuation and bifurcation software for 
delay-differential equations</i>, 2009.
</p>
<p>
<a name="man">[2]</a> R. Szalai, <i> PDDE-CONT: A continuation and bifurcation software for 
delay-differential equations</i>, 2005.
</p>
<p>
[3] R. Szalai, G. Stepan, and S. J. Hogan. 
<i>Continuation of bifurcations in periodic delay-differential equations using characteristic matrices</i>.
SIAM Journal on Scientific Computing, 28(4):1301-1317, 2006.
</p>

<p /><address><hr /><div align="center">
<table width="100%" cellspacing="0" border="0"><tr class="address">
<td width="50%">Copyright &copy; 2005, 2006, 2007 <a href="http://seis.bris.ac.uk/~rs1909">Robert Szalai</a></td>
<td width="50%" align="right"><div align="right">KNUT @PACKAGE_VERSION@</div></td>
</tr></table></div></address>
</p>
</td></tr></table>
</body></html>
