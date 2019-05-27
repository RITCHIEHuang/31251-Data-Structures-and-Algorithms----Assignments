<html>

<div class="chal-description"><h1 class="chal-title">Assignment 1 - Directed Graphs</h1> <div class="amber-display"><h2 class="amber-el amber-content">Overview</h2><p class="amber-el amber-content">Your task is to implement a directed graph class, offering a reasonably effective suite of operations, including computing spanning trees, depth and breadth first traversals, and implementing iterators.</p><h2 class="amber-el amber-content">The Code</h2><p class="amber-el amber-content">You are provided with a <code class="amber-el">directed_graph.hpp</code> file, which includes all the basic definitions you will need (and it can be done with just those - I made sure). You made add extra methods, classes, structs etc., as long as they don't interfere with the operation of the tests. You may <b class="amber-el"><u class="amber-el">not</u></b> include any further classes from the standard library in any of your marked code. You have also been provided with a <code class="amber-el">main.cpp</code> for ad-hoc testing purposes. <code class="amber-el">main.cpp</code> file does not form part of the assignment, and will not be marked. You can do anything you like with it. When the "run" button is pressed, it will compile and run <code class="amber-el">main.cpp</code>. When the "mark" button is pressed, your code will be run against the tests. Note that as this is C++, if your code causes a program crash (e.g. a segfault), the testing code <b class="amber-el"><u class="amber-el">cannot</u></b> recover, and will fail to mark your work - if you get this, make sure you fix that problem first!</p><p class="amber-el amber-content">You have also been given a copy of <code class="amber-el">test.h</code>, which is the set of tests that the code will be run against. Note that this code has had certain elements removed (things that would just give part of the solution), so it's not runnable in its current state, but it may help you understand where you program is failing, and to design your own tests to be run with <code class="amber-el">main.cpp</code>.</p><p class="amber-el amber-content">Remember to read over all the code before starting.</p><p class="amber-el amber-content">You have terminal access if you so desire it.</p><h2 class="amber-el amber-content">Directed Graphs</h2><p class="amber-el amber-content">As the abstract data structure, and the possibilities for implementing it, have been covered in the lectures, I won't repeat them here. Please refer to the lecture material for the technical details in this regard. However, don't hesitate to ask questions etc. - you're welcome to inquire, I just don't want to clutter this space up!</p><h2 class="amber-el amber-content">The <code class="amber-el">directed_graph</code> Class</h2><p class="amber-el amber-content"><code class="amber-el">directed_graph</code> is probably the most complicated C++ class you will have had to implement, but it bears a great resemblance to many of the things we've already done, so just break it down into smaller, more manageable tasks. The code itself is commented to indicate the purpose of each method. Again, to avoid clutter, I won't repeat it all here, but do not hesitate to ask if anything is unclear (there's a forum specifically for the assignment).</p><h2 class="amber-el amber-content">Marking</h2><p class="amber-el amber-content">The assignment will be marked against three components: functionality, design and style.</p><p class="amber-el amber-content"><b class="amber-el">Functionality</b> will be marked exclusively by the tests, and constitutes 50% of the total mark (note that the marks for the tests add up to 50 - so they give you the percentage you will get as well).</p><p class="amber-el amber-content"><b class="amber-el">Design</b> will be marked in your Week 9 tutorial by your tutor and constitutes 35% of the total mark. It does not depend on the functionality of your code. You may be asked questions by your tutor to help them test your understanding of your code. It will be marked qualitatively against the following rubric:</p><ul class="amber-el"><li class="amber-el amber-content"><p class="amber-el amber-content"><u class="amber-el">Pass</u> The code shows basic understanding of how to employ data<br class="amber-el">structures to achieve a goal. The design should avoid unnecessary<br class="amber-el">data structures and should make reasonable use of<br class="amber-el">iteration and recursion when appropriate.</p></li><li class="amber-el amber-content"><p class="amber-el amber-content"><u class="amber-el">Credit</u> The design shows a solid understanding of data structures and<br class="amber-el">demonstrate effective use of control structures to achieve the<br class="amber-el">program’s goals.</p></li><li class="amber-el amber-content"><p class="amber-el amber-content"><u class="amber-el">Distinction</u> The design shows a high degree of understanding of how to<br class="amber-el">use data structures to achieve a goal efficiently, and demonstrate<br class="amber-el">some evidence that the design <i class="amber-el">does not use unnecessary<br class="amber-el">resources</i>. The design should be <i class="amber-el">clean</i> and <i class="amber-el">efficient</i>.</p></li><li class="amber-el amber-content"><p class="amber-el amber-content"><u class="amber-el">High Distinction</u> The design demonstrates a high degree of understanding of<br class="amber-el">data structures and how to <i class="amber-el">efficiently</i> employ them to build<br class="amber-el">algorithms that not only meet technical goals, but support<br class="amber-el">maintenance and future development.</p></li></ul><p class="amber-el amber-content"><b class="amber-el">Style</b> will also be mark in your Week 9 tutorial by your tutor. It will be marked qualitatively against the following rubric:</p><ul class="amber-el"><li class="amber-el amber-content"><p class="amber-el amber-content"><u class="amber-el">Pass</u> The code mostly uses some formatting standard and is somewhat<br class="amber-el">readable.</p></li><li class="amber-el amber-content"><p class="amber-el amber-content"><u class="amber-el">Credit</u> The code adheres well to a formatting standard and variables<br class="amber-el">are well named.</p></li><li class="amber-el amber-content"><p class="amber-el amber-content"><u class="amber-el">Distinction</u> At least as well formatted as for Credit standards, along with<br class="amber-el">sufficient inline commenting to explain the code.</p></li><li class="amber-el amber-content"><p class="amber-el amber-content"><u class="amber-el">High Distinction</u> Excellent formatting and variable naming. Excellent, judiciously<br class="amber-el">employed comments that explain the code without just<br class="amber-el">repeating the code.</p></li></ul><p class="amber-el amber-content"><b class="amber-el">Marking Schedule</b></p><p class="amber-el amber-content">All being well, assuming you submit before the deadline and attend your Week 9 tutorial (in the class you are actually enrolled in - exceptions will only be made for cases where attendance at your enrolled tutorial is impossible), we aim to return the marks for the assignment within a week of the tutorial. However, we reserve the right to delay this schedule should technical problems arise.</p><h2 class="amber-el amber-content">Submission</h2><p class="amber-el amber-content">You will submit your work with the "mark" button on Ed. No other submissions will be accepted. You are welcome to develop your code elsewhere, if that suits your workflow, <u class="amber-el">but remember it must compile and run on Ed</u>. We are using the g++ compiler set to C++17 standard for this assessment. However code that compiles with clang++ should also work with g++ (unless you're doing something weird - so check first!).</p><p class="amber-el amber-content">You may submit as many times as you like (in fact, repeated submission is an expected part of the development process).</p><h2 class="amber-el amber-content">Due Date</h2><p class="amber-el amber-content">The assignment is due before midnight on the Friday of Week 8, which should be the 10th of May (but check that date).</p><h2 class="amber-el amber-content">Plagiarism Checking</h2><p class="amber-el amber-content">All code will be checked for student misconduct and plagiarism using specialised code similarity detection software.</p></div></div>


</html>