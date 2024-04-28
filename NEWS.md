# ProjectSVR 0.2.1

- Add `split.by` param for `MapQuery()` to project a large query Seurat object by a list of subset objects in order. 

  Now the reference mapping will be accelerated using the following R codes.
  ```R
  seu.q <- ProjectSVR(seu.q, reference, split.by = "orig.ident")
  ```

# ProjectSVR 0.2.0

- Add reference mapping wrappers and nice plot functions.

# ProjectSVR 0.1.0

- Release ProjectSVR for mapping query cells onto a well-constructed reference atlas.

