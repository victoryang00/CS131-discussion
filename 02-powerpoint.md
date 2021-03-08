DFA and NFA
===========

Introduction
------------

Finite state automata (FSA) 

RE$\leftrightarrow\epsilon$-NFA$\leftrightarrow$NFA$\leftrightarrow$DFA$\leftrightarrow$Min-DFA 

---

System as a FSA

![image-20210308174109105](/Users/yiweiyang/Library/Application Support/typora-user-images/image-20210308174109105.png)

---

Logic as a MDP NFA

Before NLP, we utilize state machine to accept languages. 

![image-20210308174251452](/Users/yiweiyang/Library/Application Support/typora-user-images/image-20210308174251452.png)

---

Wor2Vec

![image-20210308174351881](/Users/yiweiyang/Library/Application Support/typora-user-images/image-20210308174351881.png)

![image-20210308174456341](/Users/yiweiyang/Library/Application Support/typora-user-images/image-20210308174456341.png)

Examples of accepting x86 Assembly
----------------------------------

https://github.com/yangminz/bcst_csapp/

![image-20210308174724172](/Users/yiweiyang/Library/Application Support/typora-user-images/image-20210308174724172.png)

---

Assembly 模拟器

![image-20210308174757035](/Users/yiweiyang/Library/Application Support/typora-user-images/image-20210308174757035.png)

Examples of accepting yml
-------------------------

We cannot do it, because of indent can only be accepted with stack
machine.

![image-20210308174838733](/Users/yiweiyang/Library/Application Support/typora-user-images/image-20210308174838733.png)



Pumpingg Lemma
==============

Pumping Lemma: For any DFA (or NFA or regular expression) that accepts
an infinite number of strings, there is some minimum length, $M,$ such
that any string with length greater than or equal to $M$ that the
machine accepts must have the form $u x v,$ where $u, x,$ and $v$ are
strings, $x$ is not empty, the length of $u x$ is $\leq M,$ and the
machine accepts all strings of the form $u x^{n} v,$ for all $n \geq 0$.

Proof of not regular
--------------------

Let $A=\left\{1^{j} z \mid z \in\{0,1\}^{*}\right.$ and $z$ contains at
most $j$ i's, for any $\left.j \geq 1\right\} .$ Prove, by the pumping
lemma, that $A$ is not regular.\
**Answer:** 

![image-20210308175046709](/Users/yiweiyang/Library/Application Support/typora-user-images/image-20210308175046709.png)