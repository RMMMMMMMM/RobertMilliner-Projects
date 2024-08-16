using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations.Schema;
using System.ComponentModel.DataAnnotations;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Infrastructure.Models
{
    public class Review
    {

        [Key]
        public int Id { get; set; }

        [ForeignKey("RentalId")]
        public int RentalId { get; set; }

        public virtual Rental? Rental { get; set; }

        [ForeignKey("UserName")]
        [Display(Name = "Owner")]
        public AppUser? Reviewer { get; set; }

        [Required]
        public string? body { get; set; }
    }
}
