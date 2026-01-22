---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.18.1
  kernelspec:
    display_name: Python (Pyodide)
    language: python
    name: python
---

```python
import lxml
from lxml import etree
```

## Creating Elements and SubElements

```python
root = etree.Element("root")
root
```

We access the tag with the tag attribute:

```python
root.tag
```

Elements are organized in an XML tree structure. We can use the `append()` method to create a child of an element:

```python
root.append(etree.Element('child1'))
```

However, the following way is considered to be most efficient, as it assigns the child element to a variable. It uses the top-level `etree.SubElement(parent, 'child_name')` method:

```python
child2 = etree.SubElement(root, 'child2')
child3 = etree.SubElement(root, 'child3')
```

We can serialise the tree to verify it's 'XMLness'

```python
etree.tostring(root)
```

This helper function is supposed to pretty-print the XML for us, although I don't fully comprehend it (what is the `pretty_print` attribute, and why the **kwargs?

```python
def prettyprint(element, **kwargs):
    xml = etree.tostring(element, pretty_print=True, **kwargs)
    print(xml.decode(), end='')
```

```python
prettyprint(root)
```

## Elements are lists


To make access to the subelements as straight-foward as possible, the Elements class mimics the behaviour of lists as close as possible:

```python
child = root[0]
child.tag
```

```python
len(root)
```

```python
root.index(root[1]) # lxml.etree only
```

```python
children = list(root)
for child in root:
    print(child.tag)
```

```python
root.insert(0, etree.Element("child0"))
start = root[:1]
end = root[-1:]
```

```python
start[0].tag
```

```python
end[0].tag
```

To **test if an element has children**, we can use `if len(element):`, as if len == 0, this test is false and the element has no roots.


One difference from lists is that assigning an element to a different position MOVES the element, instead of copying it.

We can access the element's neighbours with `.getprevious()` and `.getnext()` 


## Elements carry attributes as a dict


xml elements have attributes. These can be created in the Element factory:

```python
root = etree.Element('root', interesting='totally')
etree.tostring(root)
```

As attributes are unordered name-value pairs, a way of dealing with them is with the dictionary-like interface of Elements:

```python
root.get('interesting')
```

```python
print(root.get('hello'))
```

```python
root.set("hello", "Huhu")
root.get('hello')
```

```python
etree.tostring(root)
```

```python
sorted(root.keys())
```

```python
for name, value in root.items():
    print(f'{name} = {value}')
```

The `attrib` attribute of an element is a real dictionary, supporting the dictionary indexing syntax of python.


## Elements contain text

```python
root = etree.Element('root')
root.text = 'Text'
root.text
```

```python
etree.tostring(root)
```

Sometimes text can surround an element, such as in the example:

`<html><body>Hello<br/>World</body></html>`

For these cases, we access the text **after** the element with the `tail` property:

```python
html = etree.Element('html')
body = etree.SubElement(html, 'body')
body.text = 'TEXT'

etree.tostring(html)
```

```python
br = etree.SubElement(body, 'br')
etree.tostring(html)
```

```python
br.tail = 'TAIL'
etree.tostring(html)
```

```python
etree.tostring(br)
```

```python
etree.tostring(br, with_tail=False)
```

```python
etree.tostring(html, method='text')
```

## Using XPath to find text


Another way to extract the text from an xml document is with the `xpath()` method, which allows you to extract the text in a "list of texts" object:

```python
text = html.xpath('//text()')
text
```

The elements of these lists are intelligent, in the way that they have a special `getparent()` method to extract what node was it's parent, returning the Element object of such parent:

```python
text[0].getparent().tag
```

```python
text[1].getparent().tag
```

```python

```
